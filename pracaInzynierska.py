from PySide6.QtWidgets import *
from PySide6.QtCore import Qt
import sys, webbrowser, os
from math import *
import math
import numpy as np
from geographiclib.geodesic import Geodesic
from geographiclib.polygonarea import PolygonArea

GRS80_a = 6378137
GRS80_e2 = 0.006694380022900788
kras_a = 6378245
kras_e2 = 0.00669342162296

def load_data_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Plik {path} nie istnieje.")
    with open(path, "r", encoding="utf-8") as f:
        content = f.read().strip()
    if not content:
        raise ValueError("Plik jest pusty.")
    if "," in content:
        raise ValueError("W pliku wykryto przecinki. Dozwolony separator dziesiętny to kropka (.).")
    content = content.replace(";", " ")
    tokens = content.split()
    try:
        numbers = [float(x) for x in tokens]
    except ValueError:
        raise ValueError("Plik zawiera wartości nienumeryczne.")
    if len(numbers) == 0:
        raise ValueError("Brak danych liczbowych w pliku.")
    if len(numbers) % 3 != 0:
        raise ValueError("Liczba wartości musi być podzielna przez 3 (X Y Z lub f l h).")
    return [numbers[i:i+3] for i in range(0, len(numbers), 3)]

def show_success(widget, operation, count, save_path, preview_lines):
    text = (
        f"SUKCES\n"
        f"Operacja: {operation}\n"
        f"Liczba przetworzonych: {count}\n"
        f"Zapisano: {save_path}\n\n"
        + "\n".join(preview_lines[:20])
    )
    widget.setText(text)

def show_error(widget, message):
    widget.setText(f"BŁĄD\nOpis: {message}")

def pk1(f,a,e2):
    N = a / np.sqrt(1 - e2 * np.sin(f)**2)
    return(N)

def XYZ2flh(X, Y, Z, a, e2):
    l = np.arctan2(Y, X)
    p = np.sqrt(X**2 + Y**2)
    f = np.arctan(Z / (p * (1 - e2)))
    while True:
        N = pk1(f,a,e2)
        h = p/np.cos(f) - N
        fs = f
        f = np.arctan(Z / (p*(1 - (e2*(N / (N + h))))))
        if np.abs(fs - f) < (0.000001/206265):
            break
    return(f,l,h)

def dms(x):
    sig = ''
    if x < 0:
        sig = '-'
        x = abs(x)
    x_deg = x * 180 / pi
    d = int(x_deg)
    m = int(60 * (x_deg - d))
    s = (x_deg - d - m/60) * 3600
    return f"{sig}{d:3d}°{abs(m):02d}'{abs(s):07.5f}\""

def dms_to_dd(d, m, s):
    sign = -1 if d < 0 else 1
    return d + sign * (m / 60.0) + sign * (s / 3600.0)

def flh2XYZ(f, l, h, a, e2):
    N = pk1(f, a, e2)
    X = (N+h) * np.cos(f) * np.cos(l)
    Y = (N+h) * np.cos(f) * np.sin(l)
    Z = (N * (1 - e2) + h) * np.sin(f)
    return(X,Y,Z)

def pkpp(f, a, e2):
    M = a * (1 - e2) / (np.sqrt(1 - e2 * np.sin(f)**2))**3
    return(M)

def kivioj(f,l,A,s,a,e2):
    n = int(s/1000)
    ds = s/n
    for i in range(n):
        M = pkpp(f,a,e2)
        N = pk1(f,a,e2)
        df = ds * np.cos(A) / M
        dA = ds * np.sin(A) * np.tan(f) / N
        fm = f + df/2
        Am = A + dA/2
        Mm = pkpp(fm,a,e2)
        Nm = pk1(fm,a,e2)
        df = (ds * np.cos(Am)) / Mm
        dA = (ds * np.sin(Am) * np.tan(fm)) / Nm 
        dl = (ds * np.sin(Am)) / (Nm * np.cos(fm))
        f = f + df
        l = l + dl 
        A = A + dA
    A = A + pi
    if A > 2*pi:
        A = A - 2*pi
    return(f,l,A)

def Vincenty(fa, la, fb, lb, a, e2):
    b = a * np.sqrt(1 - e2)
    fl = 1 - b/a
    dL = lb - la
    Ua = np.arctan((1-fl) * np.tan(fa))
    Ub = np.arctan((1-fl) * np.tan(fb))
    L = dL
    while True:
        sin_s = np.sqrt((np.cos(Ub) * np.sin(L))**2 + (np.cos(Ua) * np.sin(Ub) - np.sin(Ua) * np.cos(Ub) * np.cos(L))**2)
        cos_s = (np.sin(Ua) * np.sin(Ub)) + (np.cos(Ua) * np.cos(Ub) * np.cos(L))
        sigma = np.arctan2(sin_s, cos_s)
        sin_a = (np.cos(Ua) * np.cos(Ub) * np.sin(L)) / sin_s
        cos2_a = 1 - (sin_a)**2
        cos2sm = cos_s - ((2 * np.sin(Ua) * np.sin(Ub)) / cos2_a)
        C = fl/16 * cos2_a * (4 + fl * (4 - 3 * cos2_a))
        Ls = L
        L = dL + (1 - C) * fl * sin_a * (sigma + C * sin_s * (cos2sm + C * cos_s * (-1 + 2 * (cos2sm)**2)))
        if np.abs(Ls - L) < (0.000001/206265):
            break
    u2 = ((a**2 - b**2) / b**2) * cos2_a
    A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = (u2/1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    dsigma = B * sin_s * (cos2sm + (1/4 * B) * (cos_s * (-1 + 2 * (cos2sm)**2) - (1/6 * B * cos2sm) * (-3 + 4 * (sin_s)**2) * (-3 + 4 * (cos2sm)**2)))
    s = b * A *(sigma - dsigma)
    Aab = np.arctan2((np.cos(Ub) * np.sin(L)),((np.cos(Ua) * np.sin(Ub)) - (np.sin(Ua) * np.cos(Ub) * np.cos(L))))
    Aba = np.arctan2((np.cos(Ua) * np.sin(L)),(-np.sin(Ua) * np.cos(Ub)) + (np.cos(Ua) * np.sin(Ub) * np.cos(L))) + pi
    Aab = Aab % (2 * pi)
    Aba = Aba % (2 * pi)    
    return(s, Aab, Aba)

def sigma1(fa,a,e2):
    A0=1-e2/4-3*e2**2/64-5*e2**3/256
    A2=(3/8)*(e2+e2**2/4+15*e2**3/128)
    A4=15/256*(e2**2+3*e2**3/4)
    A6=35*e2**3/3072
    sigma=a*(A0*fa-A2*sin(2*fa)+A4*sin(4*fa)-A6*sin(6*fa))
    return(sigma)

def GKfl2xy(fa,la,l0,a,e2):
    b2=a**2*(1-e2)
    e_2=(a**2-b2)/b2
    dl=la-l0
    t=tan(fa)
    n2=e_2*cos(fa)**2
    N=pk1(fa,a,e2)
    sigma=sigma1(fa,a, e2)
    x=sigma+(dl**2/2)*N*sin(fa)*cos(fa)*(1+(dl**2/12)*(cos(fa))**2*(5-t**2+9*n2+4*n2**2)+(dl**4/360)*(cos(fa))**4*(61-58*t**2+t**4+270*n2-330*n2*t**2))
    y=dl*N*cos(fa)*(1+(dl**2/6)*(cos(fa))**2*(1-t**2+n2)+(dl**4/120)*(cos(fa))**4*(5-18*t**2+t**4+14*n2-58*n2*t**2))
    return(x,y)

def GKxy2fl(x,y,l0,a,e2):
    b2=a**2*(1-e2)
    e_2=(a**2-b2)/b2
    A0=1-e2/4-3*e2**2/64-5*e2**3/256
    fl=x/(a*A0)
    while True:
        fs=fl
        sigma=sigma1(fl,a,e2)
        fl=fl+(x-sigma)/(a*A0)
        if abs( fl-fs)<(0.000001/206265):
            break
    N=pk1(fl, a, e2)
    M=pkpp(fl, a, e2)
    t=tan(fl)
    n2=e_2*(cos(fl)**2)
    f= fl - (((y**2) * t)/(2 * M * N)) * (1 - ((y**2)/(12*N**2)) * (5 + 3 * t**2 + n2 - 9 * n2 * t**2 - 4 * n2**2) + ((y**4)/(360 * N**4)) * (61 + 90 * t**2 + 45 * t**4))
    l=l0+(y/(N*cos(fl)))*(1-((y**2)/(6*N**2))*(1+2*t**2+n2)+((y**4)/(120*N**4))*(5+28*t**2+24*t**4+6*n2+8*n2*t**2))
    return(f,l)

def f1(x, a, e2):
    A0=1-e2/4-3*e2**2/64-5*e2**3/256
    fl=x/(a*A0)
    while True:
        fs=fl
        sigma=sigma1(fl,a,e2)
        fl=fl+(x-sigma)/(a*A0)
        if abs( fl-fs)<(0.000001/206265):
            break
    return(fl)

def GK21992(x,y):
    m1992=0.9993
    x1992=x*m1992-5300000
    y1992=y*m1992+500000
    return(x1992,y1992)

def PL19922GK(x1992,y1992):
    m1992=0.9993
    x=(x1992+5300000)/m1992
    y=(y1992-500000)/m1992
    return(x,y)

def GK22000(x,y,nr):
    x2000=x*0.999923
    y2000=y*0.999923+((nr*1000000)+500000)
    return(x2000,y2000)

def PL20002GK(x2000,y2000):
    nr=int(y2000*1e-6)
    x=x2000/0.999923
    y=(y2000-((nr*1000000)+500000))/0.999923
    return(x,y)

def gamma(xgk, ygk, a, e2):
    fl = f1(xgk, a, e2)
    b2 = a**2 * (1 - e2)
    e_2 = (a**2 - b2) / b2
    t = np.tan(fl)
    n2 = e_2 * (np.cos(fl)**2)
    N = pk1(fl, a, e2)
    gamma_val = (ygk / N) * t * (1 - ((ygk**2 / (3 * (N**2))) * (1 + t**2 - n2 - 2 * n2**2)) + ((ygk**4 / (15 * (N**4))) * (2 + 5 * t**2 + 3 * t**4)))
    return gamma_val

def mgk(Xgk, Ygk, a, e2):
    f = f1(Xgk, a, e2)
    R = np.sqrt(pk1(f, a, e2) * pkpp(f, a, e2))
    m = 1 + Ygk**2 / (2 * R**2) + Ygk**4 / (24 * R**4)
    return(m)

def red_gk(xa,ya,xb,yb,l0,a,e2):
    xm = (xa + xb)/2
    ym = (ya + yb)/2
    fm,lm = GKxy2fl(xm,ym,l0,a,e2)
    Rm = np.sqrt(pkpp(fm,a,e2) * pk1(fm,a,e2))
    sgk = np.sqrt((xa - xb)**2 + (ya - yb)**2)
    r = sgk * (ya**2 + ya*yb + yb**2)/(6 * Rm**2)
    selip = sgk - r
    return(r,selip,sgk)

def red_odl1992(xa,ya,ha,xb,yb,hb,sab,l0,a,e2):
    xm = (xa + xb)/2
    ym = (ya + yb)/2
    xmgk, ymgk = PL19922GK(xm, ym)
    fm, lm = GKxy2fl(xmgk, ymgk, l0, a, e2)
    Rm = np.sqrt(pkpp(fm,a,e2) * pk1(fm,a,e2))
    s0 = sqrt(((sab**2)-(hb-ha)**2)/((1+(ha/Rm))*(1+(hb/Rm))))
    selip = (2 * Rm) * asin(s0/(2 * Rm))
    xagk, yagk = PL19922GK(xa, ya)
    xbgk, ybgk = PL19922GK(xb, yb)
    sgk = selip*(1+(((yagk**2)+(yagk*ybgk)+(ybgk**2))/(6*(Rm**2))))
    s1992 = sgk * 0.9993
    return(selip, sgk, s1992)

def red_odl2000(xa,ya,ha,xb,yb,hb,sab,l0,a,e2):
    xm = (xa + xb)/2
    ym = (ya + yb)/2
    xmgk, ymgk = PL20002GK(xm, ym)
    fm, lm = GKxy2fl(xmgk, ymgk, l0, a, e2)
    Rm = np.sqrt(pkpp(fm,a,e2) * pk1(fm,a,e2))
    s0 = sqrt(((sab**2)-(hb-ha)**2)/((1+(ha/Rm))*(1+(hb/Rm))))
    selip = (2 * Rm) * asin(s0/(2 * Rm))
    xagk, yagk = PL20002GK(xa, ya)
    xbgk, ybgk = PL20002GK(xb, yb)
    sgk = selip*(1+(((yagk**2)+(yagk*ybgk)+(ybgk**2))/(6*(Rm**2))))
    s2000 = sgk * 0.999923
    return(selip, sgk, s2000)

def delta(xa, ya, xb, yb, l0, a, e2):
    xm = (xa + xb) / 2
    ym = (ya + yb) / 2
    fm, lm = GKxy2fl(xm, ym, l0, a, e2)
    Rm = np.sqrt(pkpp(fm, a, e2) * pk1(fm, a, e2))
    deltaAB = ((xb - xa) * (2 * ya + yb) / (6 * Rm**2))
    deltaBA = ((xa - xb) * (2 * yb + ya) / (6 * Rm**2))
    return (deltaAB, deltaBA)

def alfa(xa,ya,xb,yb):
    dx = xb - xa
    dy = yb - ya
    alfa_ab = atan2(dy,dx)
    alfa_ba = alfa_ab + pi
    return(alfa_ab,alfa_ba)

def red_az(xa, ya, xb, yb, l0, a, e2):
    ga = gamma(xa, ya, a, e2)
    gb = gamma(xb, yb, a, e2)
    dab, dba = delta(xa, ya, xb, yb, l0, a, e2)
    alfa_ab, alfa_ba = alfa(xa, ya, xb, yb)
    Aab = alfa_ab + ga + dab
    Aba = alfa_ba + gb + dba
    if Aab < 0:
        Aab += 2 * pi
    if Aba < 0:
        Aba += 2 * pi
    return (Aab, Aba)

def transformacja(Kx, Ky, Kz, alfa, beta, gamma, X0, Y0, Z0, Xprim, Yprim, Zprim):
    wektor_prim = np.array([[Xprim], [Yprim], [Zprim]])
    T = np.array([[Kx, gamma, -beta], [-gamma, Ky, alfa], [beta, -alfa, Kz]])
    wektor_0 = np.array([[X0], [Y0], [Z0]])
    punkt_po_transformacji = wektor_prim + T @ wektor_prim + wektor_0
    X_bis = punkt_po_transformacji[0]
    Y_bis = punkt_po_transformacji[1]
    Z_bis = punkt_po_transformacji[2]
    return(X_bis, Y_bis, Z_bis)

def radians_dms(stopnie, minuty, sekundy):
    stopnie_dziesietne = stopnie + (minuty / 60) + (sekundy / 3600)
    radians = stopnie_dziesietne * pi / 180
    return(radians)

def transfer(x, y, z, alfa, beta, gamma, kx, ky, kz, x0, y0, z0):
    r1 = np.array([x, y, z])
    r0 = np.array([x0, y0, z0])
    M = np.array([[kx, gamma, -beta], [-gamma, ky, alfa], [beta, -alfa, kz]])
    r2 = r1 + M @ r1 + r0
    return(r2)

def polePowierzchni(filepath, a, e2):

    f_flat = 1.0 - math.sqrt(1.0 - float(e2))
    geod = Geodesic(float(a), f_flat)

    coords = []
    with open(filepath, 'r', encoding='utf-8') as fh:
        for line in fh:
            parts = line.strip().replace(',', '.').split()
            if len(parts) >= 2:
                try:
                    coords.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
 
    if len(coords) < 3:
        raise ValueError(
            f"Wielokąt musi mieć co najmniej 3 wierzchołki (wczytano {len(coords)})."
        )
    poly = PolygonArea(geod)
    for phi, lam in coords:
        poly.AddPoint(math.degrees(phi), math.degrees(lam))
 
    _, _perimeter, area = poly.Compute()
    return abs(area)

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Praca inżynierska")
        self.resize(500, 650)
        self.stack = QStackedWidget()
        self.menu_page = self.create_menu()
        self.stack.addWidget(self.menu_page)
        layout = QVBoxLayout()
        layout.addWidget(self.stack)
        self.setLayout(layout)

    def create_menu(self):
        page = QWidget()
        layout = QVBoxLayout()

        main_header = QLabel("MENU GŁÓWNE")
        main_header.setStyleSheet("font-size: 18px; font-weight: bold; color: white; padding: 10px;")
        main_header.setAlignment(Qt.AlignCenter)
        layout.addWidget(main_header)

        sub_header = QLabel("Wybierz elipsoidę referencyjną:")
        sub_header.setAlignment(Qt.AlignCenter)
        layout.addWidget(sub_header)

        self.radio_grs = QRadioButton("GRS80/WGS84")
        self.radio_kras = QRadioButton("Krasowskiego")
        self.radio_custom = QRadioButton("Własna")
        self.radio_grs.setChecked(True)

        self.radio_grs.toggled.connect(self._update_custom_visibility)
        self.radio_kras.toggled.connect(self._update_custom_visibility)
        self.radio_custom.toggled.connect(self._update_custom_visibility)

        radio_layout = QHBoxLayout()
        radio_layout.addStretch()
        radio_layout.addWidget(self.radio_grs)
        radio_layout.addWidget(self.radio_kras)
        radio_layout.addWidget(self.radio_custom)
        radio_layout.addStretch()
        layout.addLayout(radio_layout)

        self.custom_ellipsoid_box = QGroupBox("Parametry własnej elipsoidy")
        custom_grid = QGridLayout()
        custom_grid.addWidget(QLabel("a [m]:"), 0, 0)
        self.custom_a_input = QLineEdit("6378137.0")
        self.custom_a_input.setPlaceholderText("np. 6378137.0")
        custom_grid.addWidget(self.custom_a_input, 0, 1)
        custom_grid.addWidget(QLabel("e²:"), 1, 0)
        self.custom_e2_input = QLineEdit("0.00669438")
        self.custom_e2_input.setPlaceholderText("np. 0.00669438")
        custom_grid.addWidget(self.custom_e2_input, 1, 1)


        self.custom_ellipsoid_box.setLayout(custom_grid)
        self.custom_ellipsoid_box.setVisible(False)
        layout.addWidget(self.custom_ellipsoid_box)

        layout.addSpacing(10)

        task_label = QLabel("WYBIERZ ZADANIE:")
        task_label.setStyleSheet("font-weight: bold;")
        task_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(task_label)

        buttons = [
            ("Przeliczenie współrzędnych (XYZ <-> flh)", self.open_converter),
            ("Zadanie wprost", self.open_direct),
            ("Zadanie wstecz", self.open_inverse),
            ("Układy płaskie (PL-1992/PL-2000)", self.open_flat),
            ("Pole powierzchni na elipsoidzie", self.open_area),
            ("Transformacja współrzędnych przestrzennych", self.open_transformation),
            ("Pokaż punkt w Google Maps", self.open_google_maps)
        ]
        for text, func in buttons:
            btn = QPushButton(text)
            btn.setMinimumHeight(35)
            btn.clicked.connect(func)
            layout.addWidget(btn)

        layout.addStretch()
        page.setLayout(layout)
        return page

    def _update_custom_visibility(self):
        self.custom_ellipsoid_box.setVisible(self.radio_custom.isChecked())

    def get_selected_ellipsoid(self):
        if self.radio_grs.isChecked():
            return GRS80_a, GRS80_e2
        elif self.radio_kras.isChecked():
            return kras_a, kras_e2
        else:
            try:
                a_val = float(self.custom_a_input.text().replace(',', '.'))
                e2_val = float(self.custom_e2_input.text().replace(',', '.'))
            except ValueError:
                QMessageBox.critical(self, "Błąd parametrów",
                    "Niepoprawne wartości parametrów własnej elipsoidy.\n"
                    "Podaj liczby w formacie dziesiętnym (separator: kropka).")
                return None, None
            if a_val <= 0:
                QMessageBox.critical(self, "Błąd parametrów",
                    "Półoś wielka a musi być większa od 0.")
                return None, None
            if not (0 < e2_val < 1):
                QMessageBox.critical(self, "Błąd parametrów",
                    "Mimośród kwadrat e² musi należeć do przedziału (0, 1).")
                return None, None
            return a_val, e2_val

    def get_ellipsoid_name(self):
        if self.radio_grs.isChecked():
            return "GRS80/WGS84"
        elif self.radio_kras.isChecked():
            return "Krasowskiego"
        else:
            try:
                a_val = float(self.custom_a_input.text().replace(',', '.'))
                e2_val = float(self.custom_e2_input.text().replace(',', '.'))
                return f"Własna (a={a_val:.3f}, e²={e2_val:.8f})"
            except ValueError:
                return "Własna (błędne parametry)"

    def go_back(self):
        self.stack.setCurrentIndex(0)

    def _open_task(self, index, widget_class):
        a, e2 = self.get_selected_ellipsoid()
        if a is None:
            return 
        name = self.get_ellipsoid_name()
        widget = widget_class(a, e2, name, self.go_back)
        self.stack.insertWidget(index, widget)
        self.stack.setCurrentIndex(index)

    def open_converter(self):
        self._open_task(1, GeoConverter)

    def open_direct(self):
        self._open_task(2, DirectTask)

    def open_inverse(self):
        self._open_task(3, InverseTask)

    def open_flat(self):
        self._open_task(4, FlatSystems)

    def open_transformation(self):
        self._open_task(5, TransformationTask)

    def open_area(self):
        self._open_task(6, AreaTask)

    def open_google_maps(self):
        self._open_task(7, GoogleMapsTask)


class GeoConverter(QWidget):
    def __init__(self, a, e2, ell_name, go_back_callback):
        super().__init__()
        self.a = a
        self.e2 = e2
        self.go_back_callback = go_back_callback
        layout = QVBoxLayout()
        conv_header = QLabel("PRZELICZENIE WSPÓŁRZĘDNYCH")
        conv_header.setStyleSheet("font-size: 16px; font-weight: bold; color: white; padding: 5px;")
        conv_header.setAlignment(Qt.AlignCenter)
        layout.addWidget(conv_header)
        ell_info = QLabel(f"Aktywna elipsoida: {ell_name}")
        ell_info.setStyleSheet("font-weight: bold; color: white;")
        ell_info.setAlignment(Qt.AlignCenter)
        layout.addWidget(ell_info)
        layout.addSpacing(15)
        type_label = QLabel("Typ operacji:")
        type_label.setStyleSheet("color: white;")
        layout.addWidget(type_label)
        self.radio_xyz2flh = QRadioButton("XYZ -> flh")
        self.radio_flh2xyz = QRadioButton("flh -> XYZ")
        self.radio_xyz2flh.setChecked(True)
        self.conversion_group = QButtonGroup(self)
        self.conversion_group.addButton(self.radio_xyz2flh)
        self.conversion_group.addButton(self.radio_flh2xyz)
        layout.addWidget(self.radio_xyz2flh)
        layout.addWidget(self.radio_flh2xyz)
        file_layout = QHBoxLayout()
        self.file_input = QLineEdit()
        self.browse_button = QPushButton("Wybierz plik")
        self.browse_button.clicked.connect(self.choose_file)
        file_layout.addWidget(self.file_input)
        file_layout.addWidget(self.browse_button)
        layout.addLayout(file_layout)
        self.btn_convert = QPushButton("Oblicz")
        self.btn_convert.clicked.connect(self.handle_conversion)
        layout.addWidget(self.btn_convert)
        self.result_area = QTextEdit()
        self.result_area.setReadOnly(True)
        layout.addWidget(self.result_area)
        back = QPushButton("Wróć do menu głównego")
        back.clicked.connect(self.go_back_callback)
        layout.addWidget(back)
        self.setLayout(layout)

    def choose_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Wybierz plik", "", "Pliki txt (*.txt)")
        if path:
            self.file_input.setText(path)


    def handle_conversion(self):
        path = self.file_input.text()
        if not path:
            show_error(self.result_area, "Nie wybrano pliku.")
            return
        try:
            input_points = load_data_txt(path)
            if not input_points:
                raise ValueError("Brak punktów do przetworzenia.")
            output_lines = []
            for p in input_points:
                if self.radio_xyz2flh.isChecked():
                    f, l, h = XYZ2flh(p[0], p[1], p[2], self.a, self.e2)
                    output_lines.append(f"{f:.7f} {l:.7f} {h:.3f}")
                else:
                    X, Y, Z = flh2XYZ(p[0], p[1], p[2], self.a, self.e2)
                    output_lines.append(f"{X:.3f} {Y:.3f} {Z:.3f}")
            if not output_lines:
                raise ValueError("Nie wygenerowano wyników.")
            save_path = os.path.join(os.path.dirname(path), "wynikiPrzeliczenieWspolrzednych.txt")
            with open(save_path, "w", encoding="utf-8") as f_out:
                f_out.write("\n".join(output_lines))
            operation = "XYZ -> flh" if self.radio_xyz2flh.isChecked() else "flh -> XYZ"
            show_success(self.result_area, operation, len(input_points), save_path, output_lines)
        except Exception as e:
            show_error(self.result_area, str(e))

class DirectTask(QWidget):
    def __init__(self, a, e2, ell_name, go_back_callback):
        super().__init__()
        self.a = a
        self.e2 = e2
        self.go_back_callback = go_back_callback
        layout = QVBoxLayout()

        header = QLabel("ZADANIE WPROST")
        header.setStyleSheet("font-size: 16px; font-weight: bold; color: white; padding: 5px;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)

        ell_info = QLabel(f"Aktywna elipsoida: {ell_name}")
        ell_info.setStyleSheet("font-weight: bold; color: white;")
        ell_info.setAlignment(Qt.AlignCenter)
        layout.addWidget(ell_info)
        layout.addSpacing(15)

        layout.addWidget(QLabel("Typ danych w pliku:"))
        self.radio_xyz = QRadioButton("XYZ")
        self.radio_flh = QRadioButton("flh")
        self.radio_flh.setChecked(True)
        self.input_group = QButtonGroup(self)
        self.input_group.addButton(self.radio_xyz)
        self.input_group.addButton(self.radio_flh)
        layout.addWidget(self.radio_xyz)
        layout.addWidget(self.radio_flh)
        layout.addSpacing(10)

        mode_box = QGroupBox("Tryb wprowadzania azymutu i odległości")
        mode_layout = QHBoxLayout()
        self.radio_mode_gui  = QRadioButton("Ręcznie")
        self.radio_mode_file = QRadioButton("Z pliku zbiorczego")
        self.radio_mode_gui.setChecked(True)
        self.mode_group = QButtonGroup(self)
        self.mode_group.addButton(self.radio_mode_gui)
        self.mode_group.addButton(self.radio_mode_file)
        mode_layout.addWidget(self.radio_mode_gui)
        mode_layout.addWidget(self.radio_mode_file)
        mode_box.setLayout(mode_layout)
        layout.addWidget(mode_box)

        self.panel_gui = QWidget()
        gui_inner = QVBoxLayout()
        gui_inner.setContentsMargins(0, 0, 0, 0)

        az_unit_box = QGroupBox("Jednostka azymutu A1-2")
        az_unit_layout = QHBoxLayout()
        self.radio_az_grad = QRadioButton("Grady [g]")
        self.radio_az_dms  = QRadioButton("Stopnie, minuty, sekundy [° ′ ″]")
        self.radio_az_dec  = QRadioButton("Stopnie dziesiętne [°]")
        self.radio_az_grad.setChecked(True)
        self.az_unit_group = QButtonGroup(self)
        self.az_unit_group.addButton(self.radio_az_grad)
        self.az_unit_group.addButton(self.radio_az_dms)
        self.az_unit_group.addButton(self.radio_az_dec)
        az_unit_layout.addWidget(self.radio_az_grad)
        az_unit_layout.addWidget(self.radio_az_dms)
        az_unit_layout.addWidget(self.radio_az_dec)
        az_unit_box.setLayout(az_unit_layout)
        gui_inner.addWidget(az_unit_box)

        self.panel_grad = QWidget()
        grad_layout = QGridLayout()
        grad_layout.setContentsMargins(0, 4, 0, 4)
        grad_layout.addWidget(QLabel("Azymut [g]:"), 0, 0)
        self.input_azimuth = QLineEdit("0.0")
        self.input_azimuth.setPlaceholderText("np. 63.8492")
        grad_layout.addWidget(self.input_azimuth, 0, 1)
        self.panel_grad.setLayout(grad_layout)
        gui_inner.addWidget(self.panel_grad)

        self.panel_dms = QWidget()
        dms_layout = QGridLayout()
        dms_layout.setContentsMargins(0, 4, 0, 4)
        dms_layout.addWidget(QLabel("Stopnie [°]:"), 0, 0)
        self.input_az_deg = QLineEdit("0")
        self.input_az_deg.setPlaceholderText("0–359")
        dms_layout.addWidget(self.input_az_deg, 0, 1)
        dms_layout.addWidget(QLabel("Minuty [']:"), 0, 2)
        self.input_az_min = QLineEdit("0")
        self.input_az_min.setPlaceholderText("0–59")
        dms_layout.addWidget(self.input_az_min, 0, 3)
        dms_layout.addWidget(QLabel("Sekundy [\"]:"), 0, 4)
        self.input_az_sec = QLineEdit("0.0")
        self.input_az_sec.setPlaceholderText("0.0–59.999")
        dms_layout.addWidget(self.input_az_sec, 0, 5)
        self.panel_dms.setLayout(dms_layout)
        self.panel_dms.setVisible(False)
        gui_inner.addWidget(self.panel_dms)

        self.panel_dec = QWidget()
        dec_layout = QGridLayout()
        dec_layout.setContentsMargins(0, 4, 0, 4)
        dec_layout.addWidget(QLabel("Azymut [°]:"), 0, 0)
        self.input_az_dec = QLineEdit("0.0")
        self.input_az_dec.setPlaceholderText("np. 57.4642")
        dec_layout.addWidget(self.input_az_dec, 0, 1)
        self.panel_dec.setLayout(dec_layout)
        self.panel_dec.setVisible(False)
        gui_inner.addWidget(self.panel_dec)

        self.az_preview_label = QLabel("")
        self.az_preview_label.setStyleSheet("color: #aaaaaa; font-size: 10px;")
        gui_inner.addWidget(self.az_preview_label)

        dist_layout = QGridLayout()
        dist_layout.addWidget(QLabel("Odległość s [m]:"), 0, 0)
        self.input_dist = QLineEdit("1000.0")
        dist_layout.addWidget(self.input_dist, 0, 1)
        gui_inner.addLayout(dist_layout)

        self.panel_gui.setLayout(gui_inner)
        layout.addWidget(self.panel_gui)

        self.panel_batch = QWidget()
        batch_inner = QVBoxLayout()
        batch_inner.setContentsMargins(0, 0, 0, 0)

        batch_file_layout = QHBoxLayout()
        self.batch_file_input = QLineEdit()
        self.batch_file_input.setPlaceholderText("")
        self.btn_browse_batch = QPushButton("Wybierz plik")
        self.btn_browse_batch.clicked.connect(self.choose_batch_file)
        batch_file_layout.addWidget(self.batch_file_input)
        batch_file_layout.addWidget(self.btn_browse_batch)
        batch_inner.addLayout(batch_file_layout)

        self.panel_batch.setLayout(batch_inner)
        self.panel_batch.setVisible(False)
        layout.addWidget(self.panel_batch)


        self.coords_file_label = QLabel("Plik ze współrzędnymi punktów:")
        layout.addWidget(self.coords_file_label)
        file_layout = QHBoxLayout()
        self.file_input = QLineEdit()
        self.btn_browse = QPushButton("Wybierz plik")
        self.btn_browse.clicked.connect(self.choose_file)
        file_layout.addWidget(self.file_input)
        file_layout.addWidget(self.btn_browse)
        layout.addLayout(file_layout)

        self.btn_calc = QPushButton("Oblicz")
        self.btn_calc.clicked.connect(self.handle_calculation)
        layout.addWidget(self.btn_calc)

        self.result_area = QTextEdit()
        self.result_area.setReadOnly(True)
        layout.addWidget(self.result_area)

        back = QPushButton("Wróć do menu głównego")
        back.clicked.connect(self.go_back_callback)
        layout.addWidget(back)

        self.setLayout(layout)

        self.radio_az_grad.toggled.connect(self._toggle_az_unit)
        self.radio_az_dec.toggled.connect(self._toggle_az_unit)
        self.input_azimuth.textChanged.connect(self._update_az_preview)
        self.input_az_deg.textChanged.connect(self._update_az_preview)
        self.input_az_min.textChanged.connect(self._update_az_preview)
        self.input_az_sec.textChanged.connect(self._update_az_preview)
        self.input_az_dec.textChanged.connect(self._update_az_preview)
        self.radio_mode_gui.toggled.connect(self._toggle_input_mode)

    def _toggle_input_mode(self):
        gui_mode = self.radio_mode_gui.isChecked()
        self.panel_gui.setVisible(gui_mode)
        self.panel_batch.setVisible(not gui_mode)
        self.coords_file_label.setVisible(gui_mode)
        self.file_input.setVisible(gui_mode)
        self.btn_browse.setVisible(gui_mode)

    def _toggle_az_unit(self):
        self.panel_grad.setVisible(self.radio_az_grad.isChecked())
        self.panel_dms.setVisible(self.radio_az_dms.isChecked())
        self.panel_dec.setVisible(self.radio_az_dec.isChecked())
        self._update_az_preview()

    def _parse_azimuth_rad(self):
        if self.radio_az_grad.isChecked():
            try:
                val = float(self.input_azimuth.text().replace(",", "."))
            except ValueError:
                raise ValueError("Azymut w gradach: podaj liczbę dziesiętną.")
            if not (0.0 <= val < 400.0):
                raise ValueError(f"Azymut w gradach musi należeć do <0; 400). Podano: {val}")
            return val * pi / 200.0

        elif self.radio_az_dec.isChecked():
            try:
                val = float(self.input_az_dec.text().replace(",", "."))
            except ValueError:
                raise ValueError("Azymut w stopniach dziesiętnych: podaj liczbę dziesiętną.")
            if not (0.0 <= val < 360.0):
                raise ValueError(f"Azymut w stopniach musi należeć do <0; 360). Podano: {val}")
            return val * pi / 180.0

        else:
            try:
                d = float(self.input_az_deg.text().replace(",", "."))
                m = float(self.input_az_min.text().replace(",", "."))
                s = float(self.input_az_sec.text().replace(",", "."))
            except ValueError:
                raise ValueError("Azymut DMS: pola stopni, minut i sekund muszą być liczbami.")
            if not (0 <= d < 360):
                raise ValueError(f"Stopnie azymutu muszą należeć do <0; 360). Podano: {d}")
            if not (0 <= m < 60):
                raise ValueError(f"Minuty muszą należeć do <0; 60). Podano: {m}")
            if not (0.0 <= s < 60.0):
                raise ValueError(f"Sekundy muszą należeć do <0; 60). Podano: {s}")
            sign = -1 if d < 0 else 1
            deg_decimal = abs(d) + m / 60.0 + s / 3600.0
            return sign * deg_decimal * pi / 180.0

    def _update_az_preview(self):
        try:
            A_rad    = self._parse_azimuth_rad()
            deg_dec  = A_rad * 180.0 / pi
            grad_val = A_rad * 200.0 / pi
            d = int(deg_dec)
            m = int((deg_dec - d) * 60)
            s = (deg_dec - d - m / 60.0) * 3600.0

            if self.radio_az_grad.isChecked():
                self.az_preview_label.setText(
                    f"= {deg_dec:.6f}°  |  {d}° {m}′ {s:.4f}″  |  {A_rad:.7f} rad"
                )
            elif self.radio_az_dec.isChecked():
                self.az_preview_label.setText(
                    f"= {grad_val:.5f} g  |  {d}° {m}′ {s:.4f}″  |  {A_rad:.7f} rad"
                )
            else:
                self.az_preview_label.setText(
                    f"= {deg_dec:.6f}°  |  {grad_val:.5f} g  |  {A_rad:.7f} rad"
                )
        except Exception:
            self.az_preview_label.setText("")

    def choose_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Wybierz plik", "", "Pliki txt (*.txt)")
        if path:
            self.file_input.setText(path)


    def choose_batch_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Wybierz plik zbiorczy", "", "Pliki txt (*.txt)"
        )
        if path:
            self.batch_file_input.setText(path)

    def _load_batch_file(self, path):
        records = []
        with open(path, "r", encoding="utf-8") as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.replace(",", ".").split()
                if len(parts) != 5:
                    raise ValueError(
                        f"Linia {lineno}: oczekiwano 5 wartości "
                        f"(c1 c2 c3 azymut_g odleglosc_m), "
                        f"znaleziono {len(parts)}."
                    )
                try:
                    c1, c2, c3 = float(parts[0]), float(parts[1]), float(parts[2])
                    az_g       = float(parts[3])
                    s_m        = float(parts[4])
                except ValueError:
                    raise ValueError(f"Linia {lineno}: nie można przekonwertować wartości na liczby.")

                if not (0.0 <= az_g < 400.0):
                    raise ValueError(
                        f"Linia {lineno}: azymut {az_g} g poza zakresem <0; 400)."
                    )
                if s_m <= 0:
                    raise ValueError(
                        f"Linia {lineno}: odległość {s_m} m musi być > 0."
                    )

                A_rad = az_g * pi / 200.0
                records.append((c1, c2, c3, A_rad, s_m))

        if not records:
            raise ValueError("Plik zbiorczy nie zawiera żadnych danych.")
        return records

    def handle_calculation(self):
        try:
            output_lines = []

            if self.radio_mode_file.isChecked():
                batch_path = self.batch_file_input.text().strip()
                if not batch_path:
                    raise ValueError("Nie wybrano pliku zbiorczego.")

                records = self._load_batch_file(batch_path)

                for c1, c2, c3, A_rad, s_m in records:
                    if self.radio_xyz.isChecked():
                        f_s, l_s, h_s = XYZ2flh(c1, c2, c3, self.a, self.e2)
                    else:
                        f_s, l_s, h_s = c1, c2, c3
                    f_end, l_end, A_back = kivioj(f_s, l_s, A_rad, s_m, self.a, self.e2)
                    output_lines.append(
                        f"{f_end:.12f} {l_end:.12f} {A_back * 200 / pi:.5f}"
                    )

                save_path = os.path.join(os.path.dirname(batch_path), "wynikiZadanieWprost.txt")
                operation = (
                    f"Kivioj (zadanie wprost) – plik zbiorczy, "
                    f"dane: {'XYZ' if self.radio_xyz.isChecked() else 'flh'}"
                )
                n = len(records)

            else:
                path = self.file_input.text().strip()
                if not path:
                    raise ValueError("Nie wybrano pliku ze współrzędnymi.")

                A_rad = self._parse_azimuth_rad()

                try:
                    s_m = float(self.input_dist.text().replace(",", "."))
                except ValueError:
                    raise ValueError("Odległość: podaj liczbę dziesiętną.")
                if s_m <= 0:
                    raise ValueError("Odległość musi być większa od 0.")

                input_points = load_data_txt(path)
                if not input_points:
                    raise ValueError("Brak punktów do przetworzenia.")

                for p in input_points:
                    if self.radio_xyz.isChecked():
                        f_s, l_s, h_s = XYZ2flh(p[0], p[1], p[2], self.a, self.e2)
                    else:
                        f_s, l_s, h_s = p[0], p[1], p[2]
                    f_end, l_end, A_back = kivioj(f_s, l_s, A_rad, s_m, self.a, self.e2)
                    output_lines.append(
                        f"{f_end:.12f} {l_end:.12f} {A_back * 200 / pi:.5f}"
                    )

                save_path = os.path.join(os.path.dirname(path), "wynikiZadanieWprost.txt")
                if self.radio_az_grad.isChecked():
                    az_unit = "grad"
                elif self.radio_az_dec.isChecked():
                    az_unit = "stopnie dziesiętne"
                else:
                    az_unit = "DMS"
                operation = (
                    f"Kivioj (zadanie wprost), "
                    f"dane: {'XYZ' if self.radio_xyz.isChecked() else 'flh'}, "
                    f"azymut: {az_unit}"
                )
                n = len(input_points)

            if not output_lines:
                raise ValueError("Nie wygenerowano wyników.")

            with open(save_path, "w", encoding="utf-8") as f_out:
                f_out.write("\n".join(output_lines))

            show_success(self.result_area, operation, n, save_path, output_lines)

        except Exception as e:
            show_error(self.result_area, str(e))

class InverseTask(QWidget):
    def __init__(self, a, e2, ell_name, go_back_callback):
        super().__init__()
        self.a = a
        self.e2 = e2
        self.go_back_callback = go_back_callback
        layout = QVBoxLayout()
        header = QLabel("ZADANIE WSTECZ")
        header.setStyleSheet("font-size: 16px; font-weight: bold; color: white; padding: 5px;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        ell_info = QLabel(f"Aktywna elipsoida: {ell_name}")
        ell_info.setStyleSheet("font-weight: bold; color: white;")
        ell_info.setAlignment(Qt.AlignCenter)
        layout.addWidget(ell_info)
        layout.addSpacing(15)
        layout.addWidget(QLabel("Format danych wejściowych w pliku:"))
        self.radio_in_xyz = QRadioButton("XYZ")
        self.radio_in_flh = QRadioButton("flh")
        self.radio_in_flh.setChecked(True)
        self.group_input = QButtonGroup(self)
        self.group_input.addButton(self.radio_in_xyz)
        self.group_input.addButton(self.radio_in_flh)
        layout.addWidget(self.radio_in_xyz)
        layout.addWidget(self.radio_in_flh)
        layout.addSpacing(10)
        self.f_p = QLineEdit()
        btn_f = QPushButton("Wybierz plik")
        btn_f.clicked.connect(self.get_file)
        h_l = QHBoxLayout(); h_l.addWidget(self.f_p); h_l.addWidget(btn_f)
        layout.addLayout(h_l)
        btn_run = QPushButton("Oblicz")
        btn_run.clicked.connect(self.handle_inverse)
        layout.addWidget(btn_run)
        self.res = QTextEdit()
        self.res.setReadOnly(True)
        layout.addWidget(self.res)
        back = QPushButton("Wróć do menu głównego")
        back.clicked.connect(self.go_back_callback)
        layout.addWidget(back)
        self.setLayout(layout)

    def get_file(self):
        p, _ = QFileDialog.getOpenFileName(self, "Wybierz plik", "", "TXT (*.txt)")
        if p:
            self.f_p.setText(p)

    def handle_inverse(self):
        path = self.f_p.text()
        if not path:
            show_error(self.res, "Nie wybrano pliku!")
            return
        try:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Plik {path} nie istnieje.")
            output_lines = []
            processed_count = 0
            with open(path, "r", encoding="utf-8") as f:
                for line_idx, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    line_clean = line.replace(',', '.').replace(';', ' ')
                    tokens = line_clean.split()
                    if len(tokens) != 6:
                        raise ValueError(
                            f"Błąd w linii {line_idx}: Wykryto {len(tokens)} wartości.\n"
                            f"Wymagane dokładnie 6 (para punktów P1 i P2).\n"
                            f"Treść linii: {line[:30]}..."
                        )
                    try:
                        p = [float(x) for x in tokens]
                    except ValueError:
                        raise ValueError(f"Błąd w linii {line_idx}: Zawiera znaki nienumeryczne.")
                    if self.radio_in_xyz.isChecked():
                        f1, l1, h1 = XYZ2flh(p[0], p[1], p[2], self.a, self.e2)
                        f2, l2, h2 = XYZ2flh(p[3], p[4], p[5], self.a, self.e2)
                    else:
                        f1, l1, h1, f2, l2, h2 = p[0], p[1], p[2], p[3], p[4], p[5]
                    if abs(f1-f2) < 1e-11 and abs(l1-l2) < 1e-11:
                        s, Aab, Aba = 0.0, 0.0, 0.0
                    else:
                        s, Aab, Aba = Vincenty(f1, l1, f2, l2, self.a, self.e2)
                    output_lines.append(f"{s:.3f} {Aab * 200 / pi:.5f} {Aba * 200 / pi:.5f}")
                    processed_count += 1
            if not output_lines:
                raise ValueError("Plik nie zawierał żadnych danych do przetworzenia.")
            save_path = os.path.join(os.path.dirname(path), "wynikiZadanieWstecz.txt")
            with open(save_path, "w", encoding="utf-8") as f_out:
                f_out.write("\n".join(output_lines))
            show_success(self.res, "Vincenty (zadanie wstecz)", processed_count, save_path, output_lines)
        except Exception as e:
            show_error(self.res, str(e))


def auto_detect_pl2000_zone(lon_rad):
    lon_deg = lon_rad * 180.0 / pi
    meridians = {5: 15.0, 6: 18.0, 7: 21.0, 8: 24.0}
    best_zone = min(meridians, key=lambda z: abs(meridians[z] - lon_deg))
    return best_zone


class ZoneSelectionDialog(QDialog):
    def __init__(self, parent, a, e2, input_system, input_points):
        super().__init__(parent)
        self.setWindowTitle("Wybór strefy PL-2000")
        self.setMinimumWidth(380)
        self.selected_zone = None

        self._a = a
        self._e2 = e2
        self._input_system = input_system
        self._input_points = input_points

        layout = QVBoxLayout()

        desc = QLabel(
            "Układ PL-2000 podzielony jest na cztery strefy odwzorowawcze.\n"
            "Każda strefa ma swój południk osiowy:\n"
            "  Strefa 5 → λ₀ = 15°\n"
            "  Strefa 6 → λ₀ = 18°\n"
            "  Strefa 7 → λ₀ = 21°\n"
            "  Strefa 8 → λ₀ = 24°"
        )
        desc.setStyleSheet("color: #cccccc; font-size: 11px;")
        layout.addWidget(desc)
        layout.addSpacing(10)

        mode_group = QGroupBox("Tryb wyboru strefy")
        mode_layout = QVBoxLayout()

        self.radio_auto = QRadioButton("Automatyczny – na podstawie długości geograficznej punktów")
        self.radio_manual = QRadioButton("Ręczny – podaj numer strefy")
        self.radio_auto.setChecked(True)
        self.radio_auto.toggled.connect(self._toggle_mode)

        mode_layout.addWidget(self.radio_auto)
        mode_layout.addWidget(self.radio_manual)
        mode_group.setLayout(mode_layout)
        layout.addWidget(mode_group)

        self.auto_box = QGroupBox("Podgląd automatycznego wykrycia")
        auto_inner = QVBoxLayout()
        self.auto_info_label = QLabel('Kliknij "Wykryj" aby obliczyc strefę z wczytanych punktów.')
        self.auto_info_label.setWordWrap(True)
        auto_inner.addWidget(self.auto_info_label)
        btn_detect = QPushButton("Wykryj strefę automatycznie")
        btn_detect.clicked.connect(self._run_auto_detect)
        auto_inner.addWidget(btn_detect)
        self.auto_box.setLayout(auto_inner)
        layout.addWidget(self.auto_box)

        self.manual_box = QGroupBox("Numer strefy")
        manual_inner = QHBoxLayout()
        manual_inner.addWidget(QLabel("Strefa PL-2000 (5–8):"))
        self.spin_zone = QSpinBox()
        self.spin_zone.setRange(5, 8)
        self.spin_zone.setValue(6)
        manual_inner.addWidget(self.spin_zone)
        manual_inner.addStretch()
        self.manual_box.setLayout(manual_inner)
        self.manual_box.setVisible(False)
        layout.addWidget(self.manual_box)

        layout.addSpacing(8)

        btn_layout = QHBoxLayout()
        btn_ok = QPushButton("Ustaw strefę")
        btn_ok.setDefault(True)
        btn_ok.clicked.connect(self._accept)
        btn_cancel = QPushButton("Anuluj")
        btn_cancel.clicked.connect(self.reject)
        btn_layout.addStretch()
        btn_layout.addWidget(btn_ok)
        btn_layout.addWidget(btn_cancel)
        layout.addLayout(btn_layout)

        self.setLayout(layout)
        self._detected_zone = None

    def _toggle_mode(self):
        auto = self.radio_auto.isChecked()
        self.auto_box.setVisible(auto)
        self.manual_box.setVisible(not auto)

    def _run_auto_detect(self):
        """Przelicza długości geograficzne wszystkich punktów i wyznacza optymalną strefę."""
        try:
            lons_deg = []
            for p in self._input_points:
                if self._input_system == "XYZ":
                    _, l, _ = XYZ2flh(p[0], p[1], p[2], self._a, self._e2)
                    lons_deg.append(l * 180.0 / pi)
                elif self._input_system == "flh":
                    lons_deg.append(p[1] * 180.0 / pi)
                elif self._input_system == "PL-1992":
                    xgk, ygk = PL19922GK(p[0], p[1])
                    _, l = GKxy2fl(xgk, ygk, 19 * pi / 180, self._a, self._e2)
                    lons_deg.append(l * 180.0 / pi)
                elif self._input_system == "PL-2000":
                    strefa_in = int(p[1] * 1e-6)
                    xgk, ygk = PL20002GK(p[0], p[1])
                    l0 = (strefa_in * 3) * pi / 180
                    _, l = GKxy2fl(xgk, ygk, l0, self._a, self._e2)
                    lons_deg.append(l * 180.0 / pi)

            if not lons_deg:
                self.auto_info_label.setText("Brak punktów do analizy.")
                return

            lon_mean = sum(lons_deg) / len(lons_deg)
            lon_min = min(lons_deg)
            lon_max = max(lons_deg)

            meridians = {5: 15.0, 6: 18.0, 7: 21.0, 8: 24.0}

            best_zone = min(meridians, key=lambda z: abs(meridians[z] - lon_mean))
            self._detected_zone = best_zone

            zones_per_point = [min(meridians, key=lambda z: abs(meridians[z] - ld)) for ld in lons_deg]
            all_same = len(set(zones_per_point)) == 1

            info_lines = [
                f"Liczba punktów: {len(lons_deg)}",
                f"Długość geogr. min: {lon_min:.4f}°",
                f"Długość geogr. max: {lon_max:.4f}°",
                f"Długość geogr. śr.: {lon_mean:.4f}°",
                f"",
                f"➜ Wykryta strefa: {best_zone}  (λ₀ = {meridians[best_zone]}°)",
            ]
            if not all_same:
                info_lines.append(
                    "⚠ Uwaga: punkty rozrzucone po kilku strefach.\n"
                    "  Wybrano strefę dla średniej długości."
                )
            self.auto_info_label.setText("\n".join(info_lines))

        except Exception as ex:
            self.auto_info_label.setText(f"Błąd podczas wykrywania:\n{ex}")
            self._detected_zone = None

    def _accept(self):
        if self.radio_auto.isChecked():
            if self._detected_zone is None:
                QMessageBox.warning(self, "Brak wykrycia",
                    "Najpierw kliknij \"Wykryj strefę automatycznie\",\n"
                    "aby wyznaczyć strefę na podstawie danych.")
                return
            self.selected_zone = self._detected_zone
        else:
            self.selected_zone = self.spin_zone.value()
        self.accept()

class FlatSystems(QWidget):
    def __init__(self, a, e2, ell_name, go_back_callback):
        super().__init__()
        self.a = a
        self.e2 = e2
        self.go_back_callback = go_back_callback
        self.pl2000_zone = 6

        outer_layout = QVBoxLayout()

        header = QLabel("UKŁADY PŁASKIE")
        header.setStyleSheet("font-size: 16px; font-weight: bold; color: white; padding: 5px;")
        header.setAlignment(Qt.AlignCenter)
        outer_layout.addWidget(header)

        ell_info = QLabel(f"Aktywna elipsoida: {ell_name}")
        ell_info.setStyleSheet("font-weight: bold; color: white;")
        ell_info.setAlignment(Qt.AlignCenter)
        outer_layout.addWidget(ell_info)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_konwersje_tab(), "Konwersje")
        self.tabs.addTab(self._build_redukcje_tab(), "Redukcje")
        self.tabs.addTab(self._build_pole_tab(), "Pole")
        outer_layout.addWidget(self.tabs)

        back = QPushButton("Wróć do menu głównego")
        back.clicked.connect(self.go_back_callback)
        outer_layout.addWidget(back)

        self.setLayout(outer_layout)

    def _build_konwersje_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()

        self.geo_systems  = ["XYZ", "flh"]
        self.flat_systems = ["PL-1992", "PL-2000"]

        layout.addWidget(QLabel("Układ współrzędnych wejściowych:"))
        self.combo_in = QComboBox()
        self.combo_in.addItems(self.geo_systems + self.flat_systems)
        self.combo_in.currentIndexChanged.connect(self.update_output_list)
        layout.addWidget(self.combo_in)

        layout.addWidget(QLabel("Układ współrzędnych wyjściowych:"))
        self.combo_out = QComboBox()
        layout.addWidget(self.combo_out)
        self.update_output_list()

        layout.addSpacing(10)

        self.fp_konw = QLineEdit()
        btn_f = QPushButton("Wybierz plik")
        btn_f.clicked.connect(lambda: self._get_file(self.fp_konw))
        h_l = QHBoxLayout()
        h_l.addWidget(self.fp_konw)
        h_l.addWidget(btn_f)
        layout.addLayout(h_l)

        btn_calc = QPushButton("Oblicz")
        btn_calc.clicked.connect(self.handle_flat)
        layout.addWidget(btn_calc)

        self.res_konw = QTextEdit()
        self.res_konw.setReadOnly(True)
        layout.addWidget(self.res_konw)

        tab.setLayout(layout)
        return tab

    def update_output_list(self):
        selected_in = self.combo_in.currentText()
        self.combo_out.clear()
        if selected_in in self.geo_systems:
            self.combo_out.addItems(self.flat_systems)
        else:
            self.combo_out.addItems(self.geo_systems)

    def _get_file(self, line_edit):
        p, _ = QFileDialog.getOpenFileName(self, "Wybierz plik", "", "TXT (*.txt)")
        if p:
            line_edit.setText(p)

    def handle_flat(self):
        path = self.fp_konw.text()
        if not path:
            show_error(self.res_konw, "Nie wybrano pliku!")
            return
        selected_in  = self.combo_in.currentText()
        selected_out = self.combo_out.currentText()

        if selected_out == "PL-2000":
            try:
                preview_points = load_data_txt(path)
            except Exception:
                preview_points = []
            zone_dlg = ZoneSelectionDialog(self, self.a, self.e2, selected_in, preview_points)
            if zone_dlg.exec():
                self.pl2000_zone = zone_dlg.selected_zone
            else:
                return

        try:
            input_points = load_data_txt(path)
            output_lines = []
            for p in input_points:
                if selected_in == "XYZ":
                    f, l, h = XYZ2flh(p[0], p[1], p[2], self.a, self.e2)
                elif selected_in == "flh":
                    f, l, h = p[0], p[1], p[2]
                elif selected_in == "PL-1992":
                    xgk, ygk = PL19922GK(p[0], p[1])
                    f, l = GKxy2fl(xgk, ygk, 19*pi/180, self.a, self.e2)
                    h = p[2]
                elif selected_in == "PL-2000":
                    strefa_in = int(p[1] * 1e-6)
                    xgk, ygk = PL20002GK(p[0], p[1])
                    l0 = (strefa_in * 3) * pi / 180
                    f, l = GKxy2fl(xgk, ygk, l0, self.a, self.e2)
                    h = p[2]

                if selected_out == "XYZ":
                    rx, ry, rz = flh2XYZ(f, l, h, self.a, self.e2)
                elif selected_out == "flh":
                    rx, ry, rz = f, l, h
                elif selected_out == "PL-1992":
                    xgk, ygk = GKfl2xy(f, l, 19*pi/180, self.a, self.e2)
                    rx, ry = GK21992(xgk, ygk)
                    rz = h
                elif selected_out == "PL-2000":
                    l0 = (self.pl2000_zone * 3) * pi / 180
                    xgk, ygk = GKfl2xy(f, l, l0, self.a, self.e2)
                    rx, ry = GK22000(xgk, ygk, self.pl2000_zone)
                    rz = h

                if selected_out == "flh":
                    output_lines.append(f"{rx:.7f} {ry:.7f} {rz:.3f}")
                else:
                    output_lines.append(f"{rx:.3f} {ry:.3f} {rz:.3f}")

            save_path = os.path.join(os.path.dirname(path), "wynikiKonwersjiPlaskie.txt")
            with open(save_path, "w", encoding="utf-8") as f_out:
                f_out.write("\n".join(output_lines))
            show_success(self.res_konw, f"{selected_in} -> {selected_out}",
                         len(input_points), save_path, output_lines)
        except Exception as e:
            show_error(self.res_konw, str(e))

    def _build_redukcje_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()

        points_group = QGroupBox("Liczba punktów")
        points_layout = QHBoxLayout()
        self.radio_one_point  = QRadioButton("Jeden punkt  (mgk, γ)")
        self.radio_two_points = QRadioButton("Dwa punkty  (mgk, γ, δ, s, A)")
        self.radio_one_point.setChecked(True)
        points_layout.addWidget(self.radio_one_point)
        points_layout.addWidget(self.radio_two_points)
        points_group.setLayout(points_layout)
        layout.addWidget(points_group)

        layout.addWidget(QLabel("Układ współrzędnych wejściowych:"))
        self.combo_red_sys = QComboBox()
        self.combo_red_sys.addItems(["PL-1992", "PL-2000"])
        layout.addWidget(self.combo_red_sys)

        self.fp_red = QLineEdit()
        btn_red = QPushButton("Wybierz plik")
        btn_red.clicked.connect(lambda: self._get_file(self.fp_red))
        h_l = QHBoxLayout()
        h_l.addWidget(self.fp_red)
        h_l.addWidget(btn_red)
        layout.addLayout(h_l)

        btn_calc = QPushButton("Oblicz")
        btn_calc.clicked.connect(self.handle_redukcje)
        layout.addWidget(btn_calc)

        self.res_red = QTextEdit()
        self.res_red.setReadOnly(True)
        layout.addWidget(self.res_red)

        tab.setLayout(layout)
        return tab

    def _load_single_points_txt(self, path):
        """Wczytuje pojedyncze punkty: każda linia to  x y  lub  x y h."""
        if not os.path.exists(path):
            raise FileNotFoundError(f"Plik {path} nie istnieje.")
        with open(path, "r", encoding="utf-8") as f:
            content = f.read().strip()
        if not content:
            raise ValueError("Plik jest pusty.")
        if "," in content:
            raise ValueError("Wykryto przecinki. Dozwolony separator dziesiętny to kropka.")
        content = content.replace(";", " ")
        points = []
        for i, line in enumerate(content.splitlines(), 1):
            line = line.strip()
            if not line:
                continue
            tokens = line.split()
            if len(tokens) < 2:
                raise ValueError(f"Linia {i}: potrzeba co najmniej 2 wartości (x y).")
            try:
                vals = [float(t) for t in tokens[:2]]
            except ValueError:
                raise ValueError(f"Linia {i}: wartości nienumeryczne.")
            points.append((vals[0], vals[1]))
        if not points:
            raise ValueError("Brak danych w pliku.")
        return points

    def _load_pairs_txt(self, path):
        """Wczytuje pary punktów: każda linia to  xa ya ha xb yb hb."""
        if not os.path.exists(path):
            raise FileNotFoundError(f"Plik {path} nie istnieje.")
        with open(path, "r", encoding="utf-8") as f:
            content = f.read().strip()
        if not content:
            raise ValueError("Plik jest pusty.")
        if "," in content:
            raise ValueError("Wykryto przecinki. Dozwolony separator dziesiętny to kropka.")
        content = content.replace(";", " ")
        pairs = []
        for i, line in enumerate(content.splitlines(), 1):
            line = line.strip()
            if not line:
                continue
            tokens = line.split()
            if len(tokens) != 6:
                raise ValueError(
                    f"Linia {i}: wykryto {len(tokens)} wartości. "
                    "Wymagane dokładnie 6: xa ya ha xb yb hb."
                )
            try:
                vals = [float(t) for t in tokens]
            except ValueError:
                raise ValueError(f"Linia {i}: wartości nienumeryczne.")
            pairs.append((vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]))
        if not pairs:
            raise ValueError("Brak danych w pliku.")
        return pairs

    def handle_redukcje(self):
        path = self.fp_red.text()
        if not path:
            show_error(self.res_red, "Nie wybrano pliku!")
            return
        sys_in = self.combo_red_sys.currentText()

        try:
            if self.radio_one_point.isChecked():
                points = self._load_single_points_txt(path)
                output_lines = []
                header_line = "  Nr            x [m]          y [m]             mgk    gamma [rad]  gamma [°'\"]"
                output_lines.append(header_line)
                output_lines.append("-" * 95)

                for idx, (x, y) in enumerate(points, 1):
                    if sys_in == "PL-1992":
                        xgk, ygk = PL19922GK(x, y)
                    else:
                        xgk, ygk = PL20002GK(x, y)

                    m_val = mgk(xgk, ygk, self.a, self.e2)
                    g_val = gamma(xgk, ygk, self.a, self.e2)
                    g_dms_str = dms(g_val)

                    output_lines.append(
                        f"{idx:>5}  {x:>14.3f}  {y:>14.3f}  "
                        f"{m_val:>14.8f}  {g_val:>13.7f}  {g_dms_str}"
                    )

                save_path = os.path.join(os.path.dirname(path), "wynikiRedukcje.txt")
                with open(save_path, "w", encoding="utf-8") as f_out:
                    f_out.write("\n".join(output_lines))
                show_success(self.res_red, "Redukcje – jeden punkt (mgk, γ)",
                             len(points), save_path, output_lines)

            else:
                pairs = self._load_pairs_txt(path)
                output_lines = []
                col_names = (
                    f"{'Para':>5}  "
                    f"{'mgk_A':>12}  {'gamma_A[rad]':>13}  {'gamma_A[dms]':>22}  "
                    f"{'mgk_B':>12}  {'gamma_B[rad]':>13}  {'gamma_B[dms]':>22}  "
                    f"{'sAB_Vinc[m]':>13}  {'s_elip[m]':>12}  {'s_gk[m]':>12}  {'s_uklad[m]':>12}  "
                    f"{'deltaAB[rad]':>13}  {'deltaBA[rad]':>13}  "
                    f"{'alfa_AB[rad]':>13}  {'alfa_BA[rad]':>13}  "
                    f"{'Aab_red[rad]':>13}  {'Aab_red[g]':>12}  "
                    f"{'Aba_red[rad]':>13}  {'Aba_red[g]':>12}"
                )
                output_lines.append(col_names)
                output_lines.append("-" * 245)

                for idx, (xa, ya, hA, xb, yb, hB) in enumerate(pairs, 1):
                    if sys_in == "PL-1992":
                        xagk, yagk = PL19922GK(xa, ya)
                        xbgk, ybgk = PL19922GK(xb, yb)
                        l0 = 19 * pi / 180

                        fA, lA = GKxy2fl(xagk, yagk, l0, self.a, self.e2)
                        fB, lB = GKxy2fl(xbgk, ybgk, l0, self.a, self.e2)
 
                        sAB_vinc, _, _ = Vincenty(fA, lA, fB, lB, self.a, self.e2)

                        selip, sgk, s_uklad = red_odl1992(
                            xa, ya, hA, xb, yb, hB, sAB_vinc, l0, self.a, self.e2
                        )
                    else:
                        strefa = int(ya * 1e-6)
                        l0 = strefa * 3 * pi / 180
                        xagk, yagk = PL20002GK(xa, ya)
                        xbgk, ybgk = PL20002GK(xb, yb)

                        fA, lA = GKxy2fl(xagk, yagk, l0, self.a, self.e2)
                        fB, lB = GKxy2fl(xbgk, ybgk, l0, self.a, self.e2)

                        sAB_vinc, _, _ = Vincenty(fA, lA, fB, lB, self.a, self.e2)

                        selip, sgk, s_uklad = red_odl2000(
                            xa, ya, hA, xb, yb, hB, sAB_vinc, l0, self.a, self.e2
                        )

                    mA = mgk(xagk, yagk, self.a, self.e2)
                    mB = mgk(xbgk, ybgk, self.a, self.e2)

                    gA     = gamma(xagk, yagk, self.a, self.e2)
                    gB     = gamma(xbgk, ybgk, self.a, self.e2)
                    gA_str = dms(gA)
                    gB_str = dms(gB)

                    dAB, dBA = delta(xagk, yagk, xbgk, ybgk, l0, self.a, self.e2)
                    alfa_ab, alfa_ba = alfa(xagk, yagk, xbgk, ybgk)

                    Aab_red, Aba_red = red_az(xagk, yagk, xbgk, ybgk, l0, self.a, self.e2)

                    Aab_red = Aab_red % (2 * pi)
                    Aba_red = Aba_red % (2 * pi)
                    Aab_red_g = Aab_red * 200 / pi
                    Aba_red_g = Aba_red * 200 / pi

                    output_lines.append(
                        f"{idx:>5}  "
                        f"{mA:>12.8f}  {gA:>13.7f}  {gA_str:>22}  "
                        f"{mB:>12.8f}  {gB:>13.7f}  {gB_str:>22}  "
                        f"{sAB_vinc:>13.4f}  {selip:>12.4f}  {sgk:>12.4f}  {s_uklad:>12.4f}  "
                        f"{dAB:>13.7f}  {dBA:>13.7f}  "
                        f"{alfa_ab:>13.7f}  {alfa_ba:>13.7f}  "
                        f"{Aab_red:>13.7f}  {Aab_red_g:>12.5f}  "
                        f"{Aba_red:>13.7f}  {Aba_red_g:>12.5f}"
                    )

                save_path = os.path.join(os.path.dirname(path), "wynikiRedukcje.txt")
                with open(save_path, "w", encoding="utf-8") as f_out:
                    f_out.write("\n".join(output_lines))
                show_success(self.res_red,
                             "Redukcje – dwa punkty (mgk, γ, redukcja długości i azymutu, δ)",
                             len(pairs), save_path, output_lines)

        except Exception as e:
            show_error(self.res_red, str(e))


    def _build_pole_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Układ współrzędnych wejściowych:"))
        self.combo_pole_sys = QComboBox()
        self.combo_pole_sys.addItems(["PL-1992", "PL-2000"])
        layout.addWidget(self.combo_pole_sys)

        self.fp_pole = QLineEdit()
        btn_pole = QPushButton("Wybierz plik")
        btn_pole.clicked.connect(lambda: self._get_file(self.fp_pole))
        h_l = QHBoxLayout()
        h_l.addWidget(self.fp_pole)
        h_l.addWidget(btn_pole)
        layout.addLayout(h_l)

        btn_calc = QPushButton("Oblicz")
        btn_calc.clicked.connect(self.handle_pole)
        layout.addWidget(btn_calc)

        self.res_pole = QTextEdit()
        self.res_pole.setReadOnly(True)
        layout.addWidget(self.res_pole)

        tab.setLayout(layout)
        return tab

    def _load_points_xy_txt(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Plik {path} nie istnieje.")
        with open(path, "r", encoding="utf-8") as f:
            content = f.read().strip()
        if not content:
            raise ValueError("Plik jest pusty.")
        if "," in content:
            raise ValueError("Wykryto przecinki. Dozwolony separator dziesiętny to kropka.")
        content = content.replace(";", " ")
        points = []
        for i, line in enumerate(content.splitlines(), 1):
            line = line.strip()
            if not line:
                continue
            tokens = line.split()
            if len(tokens) == 2:
                try:
                    vals = [float(t) for t in tokens]
                except ValueError:
                    raise ValueError(f"Linia {i}: wartości nienumeryczne.")
                points.append((vals[0], vals[1]))
            elif len(tokens) == 3:
                try:
                    vals = [float(t) for t in tokens]
                except ValueError:
                    raise ValueError(f"Linia {i}: wartości nienumeryczne.")
                points.append((vals[0], vals[1]))
            else:
                raise ValueError(
                    f"Linia {i}: wykryto {len(tokens)} wartości. "
                    "Wymagane 2 (x y) lub 3 (x y h)."
                )
        if len(points) < 3:
            raise ValueError("Potrzeba co najmniej 3 punktów do obliczenia pola.")
        return points

    def _pole_gaussa(self, points):
        """Pole Gaussa ze współrzędnych płaskich."""
        n = len(points)
        area = 0.0
        for i in range(n):
            xi, yi = points[i]
            xn, yn = points[(i + 1) % n]
            area += xi * yn - xn * yi
        return abs(area) / 2.0

    def handle_pole(self):
        path = self.fp_pole.text()
        if not path:
            show_error(self.res_pole, "Nie wybrano pliku!")
            return
        sys_in = self.combo_pole_sys.currentText()
        try:
            points_in = self._load_points_xy_txt(path)

            points_gk = []
            m_gk_list = []  
            m_ukl_list = [] 

            if sys_in == "PL-1992":
                l0_sys = 19 * pi / 180
                m0 = 0.9993
            else:
                strefa = int(points_in[0][1] * 1e-6)
                l0_sys = strefa * 3 * pi / 180
                m0 = 0.999923

            for (x, y) in points_in:
                if sys_in == "PL-1992":
                    xgk, ygk = PL19922GK(x, y)
                else:
                    xgk, ygk = PL20002GK(x, y)

                f, l = GKxy2fl(xgk, ygk, l0_sys, self.a, self.e2)


                m_g = mgk(xgk, ygk, self.a, self.e2)
                m_gk_list.append(m_g)
                m_ukl_list.append(m_g * m0)
                points_gk.append((xgk, ygk))


            pole_plask = self._pole_gaussa(points_in)
            pole_gk = self._pole_gaussa(points_gk)
            
            m_sr_gk = sum(m_gk_list) / len(m_gk_list)
            m_sr_ukl = sum(m_ukl_list) / len(m_ukl_list)

            pole_elip_z_gk = pole_gk / (m_sr_gk ** 2)
            pole_elip_z_ukl = pole_plask / (m_sr_ukl ** 2)
            



            lines = []
            lines.append(f"Układ wejściowy: {sys_in}")
            lines.append(f"Liczba punktów: {len(points_in)}")
            lines.append(f"Mnożnik skali układu m0 = {m0}")
            lines.append("")
            lines.append(
                f"{'Punkt':>6}  {'m_gk':>12}  {'m_uklad':>12}  "
                f"{'Lzniek [cm/km]':>16}  {'Pzniek [m²/ha]':>16}"
            )
            lines.append("-" * 72)

            for i, (m_g, m_u) in enumerate(zip(m_gk_list, m_ukl_list), 1):
                L_zn = (m_u - 1) * 100_000
                P_zn = (m_u ** 2 - 1) * 10_000
                x_in, y_in = points_in[i - 1]
                lines.append(
                    f"{i:>6}  {m_g:>12.8f}  {m_u:>12.8f}  "
                    f"{L_zn:>16.4f}  {P_zn:>16.4f}"
                )

            m_sr_gk_val  = sum(m_gk_list)  / len(m_gk_list)
            m_sr_ukl_val = sum(m_ukl_list) / len(m_ukl_list)
            L_zn_sr = (m_sr_ukl_val - 1) * 100_000
            P_zn_sr = (m_sr_ukl_val ** 2 - 1) * 10_000

            lines.append("-" * 72)
            lines.append(
                f"{'Średnia':>6}  {m_sr_gk_val:>12.8f}  {m_sr_ukl_val:>12.8f}  "
                f"{L_zn_sr:>16.4f}  {P_zn_sr:>16.4f}"
            )
            lines.append("")
            lines.append("=" * 72)
            lines.append("POLE POWIERZCHNI")
            lines.append("=" * 72)
            lines.append(f"  Pole płaskie (Gaussa, wsp. układu) [m²]:    {pole_plask:.3f}")
            lines.append(f"  Pole płaskie (Gaussa, wsp. układu) [ha]:    {pole_plask / 10_000:.4f}")
            lines.append(f"  Pole płaskie (Gaussa, wsp. układu) [km²]:   {pole_plask / 1_000_000:.6f}")
            lines.append("")
            lines.append(f"  Pole GK (Gaussa, wsp. GK)          [m²]:    {pole_gk:.3f}")
            lines.append("")
            lines.append(f"  Pole na elipsoidzie (z wsp. GK)    [m²]:    {pole_elip_z_gk:.3f}")
            lines.append(f"  Pole na elipsoidzie (z wsp. GK)    [ha]:    {pole_elip_z_gk / 10_000:.4f}")
            lines.append(f"  Pole na elipsoidzie (z wsp. GK)    [km²]:   {pole_elip_z_gk / 1_000_000:.6f}")
            lines.append("")
            lines.append(f"  Pole na elipsoidzie (z wsp. ukł.)  [m²]:    {pole_elip_z_ukl:.3f}")
            lines.append(f"  Pole na elipsoidzie (z wsp. ukł.)  [ha]:    {pole_elip_z_ukl / 10_000:.4f}")
            lines.append(f"  Pole na elipsoidzie (z wsp. ukł.)  [km²]:   {pole_elip_z_ukl / 1_000_000:.6f}")
            lines.append("")
            lines.append(f"  Użyta średnia skala GK  m_sr_GK  = {m_sr_gk_val:.8f}")
            lines.append(f"  Użyta średnia skala ukł m_sr_ukl = {m_sr_ukl_val:.8f}")

            result_text = "\n".join(lines)
            self.res_pole.setText(result_text)

            save_path = os.path.join(os.path.dirname(path), "wynikiPoleGK.txt")
            with open(save_path, "w", encoding="utf-8") as f_out:
                f_out.write(result_text)

            show_success(self.res_pole, f"Pole i zniekształcenia ({sys_in})",
                         len(points_in), save_path, lines[:20])
        except Exception as e:
            show_error(self.res_pole, str(e))

class TransformationTask(QWidget):
    def __init__(self, a, e2, ell_name, go_back_callback):
        super().__init__()
        self.a = a
        self.e2 = e2
        self.go_back_callback = go_back_callback
        layout = QVBoxLayout()
        header = QLabel("TRANSFORMACJA")
        header.setStyleSheet("font-size: 16px; font-weight: bold; color: white; padding: 5px;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        ell_info = QLabel(f"Aktywna elipsoida: {ell_name}")
        ell_info.setStyleSheet("font-weight: bold; color: white;")
        ell_info.setAlignment(Qt.AlignCenter)
        layout.addWidget(ell_info)
        layout.addSpacing(15)
        layout.addWidget(QLabel("Źródło parametrów transformacji:"))
        self.radio_param_gui = QRadioButton("Podaj parametry")
        self.radio_param_file = QRadioButton("Pobierz parametry z pliku (Format: X' Y' Z' Kx Ky Kz alfa beta gamma X0 Y0 Z0)")
        self.radio_param_gui.setChecked(True)
        self.radio_param_gui.toggled.connect(self.toggle_params_visibility)
        layout.addWidget(self.radio_param_gui)
        layout.addWidget(self.radio_param_file)
        self.params_container = QGroupBox("Parametry transformacji")
        params_grid = QGridLayout()
        self.inputs = {}
        labels = [
            ("Kx", 0, 0), ("Ky", 0, 2), ("Kz", 0, 4),
            ("alfa [rad]", 1, 0), ("beta [rad]", 1, 2), ("gamma [rad]", 1, 4),
            ("X0 [m]", 2, 0), ("Y0 [m]", 2, 2), ("Z0 [m]", 2, 4)
        ]
        for name, row, col in labels:
            params_grid.addWidget(QLabel(name + ":"), row, col)
            edit = QLineEdit("0.0")
            params_grid.addWidget(edit, row, col + 1)
            self.inputs[name.split()[0]] = edit
        self.params_container.setLayout(params_grid)
        layout.addWidget(self.params_container)
        layout.addWidget(QLabel(""))
        file_layout = QHBoxLayout()
        self.file_input = QLineEdit()
        self.btn_browse = QPushButton("Wybierz plik")
        self.btn_browse.clicked.connect(self.choose_file)
        file_layout.addWidget(self.file_input)
        file_layout.addWidget(self.btn_browse)
        layout.addLayout(file_layout)
        self.btn_run = QPushButton("Oblicz")
        self.btn_run.clicked.connect(self.handle_transformation)
        layout.addWidget(self.btn_run)
        self.result_area = QTextEdit()
        self.result_area.setReadOnly(True)
        layout.addWidget(self.result_area)
        back = QPushButton("Wróć do menu głównego")
        back.clicked.connect(self.go_back_callback)
        layout.addWidget(back)
        self.setLayout(layout)

    def toggle_params_visibility(self):
        self.params_container.setEnabled(self.radio_param_gui.isChecked())

    def choose_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Wybierz plik", "", "Pliki txt (*.txt)")
        if path:
            self.file_input.setText(path)

    def handle_transformation(self):
        path = self.file_input.text()
        if not path:
            self.result_area.setText("Błąd: Nie wybrano pliku!")
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                content = f.read().replace(',', '.').replace(';', ' ')
            d = [float(x) for x in content.split()]
            output_lines = []
            if self.radio_param_gui.isChecked():
                kx = float(self.inputs["Kx"].text())
                ky = float(self.inputs["Ky"].text())
                kz = float(self.inputs["Kz"].text())
                alfa = float(self.inputs["alfa"].text())
                beta = float(self.inputs["beta"].text())
                gamma = float(self.inputs["gamma"].text())
                x0 = float(self.inputs["X0"].text())
                y0 = float(self.inputs["Y0"].text())
                z0 = float(self.inputs["Z0"].text())
                if len(d) % 3 != 0:
                    raise ValueError("Plik musi zawierać współrzędne punktów (wielokrotność 3 liczb).")
                for i in range(0, len(d), 3):
                    xb, yb, zb = transformacja(kx, ky, kz, alfa, beta, gamma, x0, y0, z0, d[i], d[i+1], d[i+2])
                    output_lines.append(f"{float(xb):.3f} {float(yb):.3f} {float(zb):.3f}")
            else:
                if len(d) % 12 != 0:
                    raise ValueError("Format pliku z parametrami musi mieć 12 liczb na punkt!")
                for i in range(0, len(d), 12):
                    xb, yb, zb = transformacja(d[i+3], d[i+4], d[i+5], d[i+6], d[i+7], d[i+8], d[i+9], d[i+10], d[i+11], d[i], d[i+1], d[i+2])
                    output_lines.append(f"{float(xb):.3f} {float(yb):.3f} {float(zb):.3f}")
            save_path = os.path.join(os.path.dirname(path), "wynikiTransformacja.txt")
            with open(save_path, "w", encoding="utf-8") as f_out:
                f_out.write("\n".join(output_lines))
            operation = "Transformacja współrzędnych przestrzennych"
            count = len(output_lines)
            show_success(self.result_area, operation, count, save_path, output_lines)
        except Exception as e:
            show_error(self.result_area, str(e))

class AreaTask(QWidget):
    """
    Zakładka 'Pole powierzchni na elipsoidzie'.

    Wczytuje plik z współrzędnymi (φ, λ w radianach)
    i oblicza pole geodezyjnego wielokąta metodą Karneya (2013) /
    Sjöberga (2006) z dokładnością ~15 cyfr.
    """

    def __init__(self, a, e2, ell_name, go_back_callback):
        super().__init__()
        self.a                = a
        self.e2               = e2
        self.go_back_callback = go_back_callback

        layout = QVBoxLayout()

        header = QLabel("POLE POWIERZCHNI NA ELIPSOIDZIE")
        header.setStyleSheet(
            "font-size: 16px; font-weight: bold; color: white; padding: 5px;"
        )
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)

        ell_info = QLabel(f"Aktywna elipsoida: {ell_name}")
        ell_info.setStyleSheet("font-weight: bold; color: white;")
        ell_info.setAlignment(Qt.AlignCenter)
        layout.addWidget(ell_info)


        layout.addSpacing(10)

        fmt_layout = QVBoxLayout()
        
        fmt_label = QLabel("Format danych wejściowych w pliku:")
        fmt_layout.addWidget(fmt_label)

        self.radio_flh = QRadioButton("flh")
        self.radio_flh.setChecked(True)
        fmt_layout.addWidget(self.radio_flh)
        
        layout.addLayout(fmt_layout)

        layout.addSpacing(10)

        file_layout = QHBoxLayout()
        self.file_input = QLineEdit()
        self.file_input.setPlaceholderText("")
        self.btn_browse = QPushButton("Wybierz plik")
        self.btn_browse.clicked.connect(self.choose_file)
        file_layout.addWidget(self.file_input)
        file_layout.addWidget(self.btn_browse)
        layout.addLayout(file_layout)

        layout.addSpacing(5)
        self.btn_calc = QPushButton("Oblicz")
        self.btn_calc.setMinimumHeight(36)
        self.btn_calc.clicked.connect(self.handle_calculation)
        layout.addWidget(self.btn_calc)

        self.result_area = QTextEdit()
        self.result_area.setReadOnly(True)
        self.result_area.setMinimumHeight(180)
        layout.addWidget(self.result_area)

        back = QPushButton("Wróć do menu głównego")
        back.clicked.connect(self.go_back_callback)
        layout.addWidget(back)

        self.setLayout(layout)

    def choose_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Wybierz plik", "", "Pliki txt (*.txt)")
        if path:
            self.file_input.setText(path)


    def _load_points_as_fl(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Plik nie istnieje: {path}")

        points = []

        with open(path, 'r', encoding='utf-8') as fh:
            for lineno, raw in enumerate(fh, 1):
                line = raw.strip().replace(';', ' ').replace(',', '.')
                if not line:
                    continue
                try:
                    vals = [float(x) for x in line.split()]
                except ValueError:
                    raise ValueError(
                        f"Linia {lineno}: nie można odczytać wartości liczbowych."
                    )

                if len(vals) < 2:
                    raise ValueError(
                        f"Linia {lineno}: wymagane co najmniej dwie wartości (φ λ)."
                    )
                
                phi, lam = vals[0], vals[1]
                points.append((phi, lam))

        if len(points) < 3:
            raise ValueError(
                f"Wielokąt musi mieć co najmniej 3 wierzchołki "
                f"(wczytano {len(points)})."
            )
        return points


    def handle_calculation(self):
        path = self.file_input.text().strip()
        if not path:
            show_error(
                self.result_area,
                "Nie wybrano pliku — kliknij 'Wybierz plik'."
            )
            return

        try:
            processed = self._load_points_as_fl(path)

            dir_path = os.path.dirname(path) or "."
            tmp_path = os.path.join(dir_path, "_tmp_pole_.txt")
            with open(tmp_path, 'w', encoding='utf-8') as tf:
                for phi, lam in processed:
                    tf.write(f"{phi:.15f} {lam:.15f}\n")

            area_m2 = polePowierzchni(tmp_path, self.a, self.e2)

            try:
                os.remove(tmp_path)
            except OSError:
                pass

            area_ha  = area_m2 / 10_000.0
            area_km2 = area_m2 / 1_000_000.0

            output_lines = [
                f"Elipsoida: a = {self.a} m,  e\u00b2 = {self.e2}",
                f"Liczba wierzcho\u0142k\u00f3w: {len(processed)}",
                "",
                f"Pole [m\u00b2]  :  {area_m2:.4f}",
                f"Pole [ha]  :  {area_ha:.6f}",
                f"Pole [km\u00b2] :  {area_km2:.8f}",
            ]

            save_path = os.path.join(dir_path, "wynikPolePowierzchni.txt")
            with open(save_path, 'w', encoding='utf-8') as fout:
                fout.write("\n".join(output_lines))

            show_success(
                self.result_area,
                "Pole powierzchni na elipsoidzie",
                len(processed),
                save_path,
                output_lines,
            )

        except ImportError as exc:
            show_error(
                self.result_area,
                str(exc) + "\n\nUruchom w terminalu:\n  pip install geographiclib"
            )
        except Exception as exc:
            show_error(self.result_area, str(exc))

class GoogleMapsTask(QWidget):
    def show_info(self, message):
        self.status_area.setText(f"INFO\n{message}")

    def __init__(self, a, e2, ell_name, go_back_callback):
        super().__init__()
        self.a = a
        self.e2 = e2
        self.go_back_callback = go_back_callback
        layout = QVBoxLayout()
        header = QLabel("LOKALIZACJA W GOOGLE MAPS")
        header.setStyleSheet("font-size: 16px; font-weight: bold; color: white; padding: 5px;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        ell_info = QLabel(f"Aktywna elipsoida: {ell_name}")
        ell_info.setStyleSheet("font-weight: bold; color: white;")
        ell_info.setAlignment(Qt.AlignCenter)
        layout.addWidget(ell_info)
        format_group = QGroupBox("Format danych:")
        format_layout = QVBoxLayout()
        self.radio_xyz = QRadioButton("XYZ")
        self.radio_rad = QRadioButton("flh (radiany)")
        self.radio_deg = QRadioButton("flh (stopnie dziesiętne)")
        self.radio_dms = QRadioButton("flh (° ′ ″)")
        self.radio_deg.setChecked(True)
        format_layout.addWidget(self.radio_xyz)
        format_layout.addWidget(self.radio_rad)
        format_layout.addWidget(self.radio_deg)
        format_layout.addWidget(self.radio_dms)
        format_group.setLayout(format_layout)
        layout.addWidget(format_group)
        manual_box = QGroupBox("Wprowadź współrzędne")
        manual_layout = QVBoxLayout()
        self.coords_input = QLineEdit()
        manual_layout.addWidget(self.coords_input)
        manual_box.setLayout(manual_layout)
        layout.addWidget(manual_box)
        file_box = QGroupBox("Załaduj z pliku .txt:")
        file_layout = QHBoxLayout()
        self.file_input = QLineEdit()
        btn_browse = QPushButton("Wybierz plik")
        btn_browse.clicked.connect(self.choose_file)
        file_layout.addWidget(self.file_input)
        file_layout.addWidget(btn_browse)
        file_box.setLayout(file_layout)
        layout.addWidget(file_box)
        btn_open = QPushButton("OTWÓRZ W GOOGLE MAPS")
        btn_open.setMinimumHeight(40)
        btn_open.setStyleSheet("background-color: #27AE60; color: white; font-weight: bold;")
        btn_open.clicked.connect(self.open_map)
        layout.addWidget(btn_open)
        self.status_area = QTextEdit()
        self.status_area.setReadOnly(True)
        self.status_area.setMaximumHeight(80)
        layout.addWidget(self.status_area)
        back = QPushButton("Wróć do menu głównego")
        back.clicked.connect(self.go_back_callback)
        layout.addWidget(back)
        self.setLayout(layout)

    def choose_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Wybierz plik", "", "Pliki txt (*.txt)")
        if path:
            self.file_input.setText(path)


    def open_map(self):
        try:
            path = self.file_input.text()
            manual_text = self.coords_input.text()
            multi_point_warning = False
            if path:
                pts = load_data_txt(path)
                if not pts:
                    raise ValueError("Plik jest pusty!")
                if len(pts) > 1:
                    multi_point_warning = True
                p = pts[0]
            elif manual_text:
                p = [float(x) for x in manual_text.replace(',', '.').split()]
                if len(p) < 2:
                    raise ValueError("Podaj odpowiednią liczbę wartości!")
            else:
                return self.status_area.setText("BŁĄD: Brak danych!")
            if self.radio_xyz.isChecked():
                if len(p) < 3:
                    raise ValueError("XYZ wymaga 3 współrzędnych!")
                f, l, _ = XYZ2flh(p[0], p[1], p[2], self.a, self.e2)
                phi, lam = degrees(f), degrees(l)
            elif self.radio_rad.isChecked():
                phi, lam = degrees(p[0]), degrees(p[1])
            elif self.radio_dms.isChecked():
                if len(p) < 6:
                    raise ValueError("Format DMS wymaga 6 wartości (° ′ ″ ° ′ ″)!")
                phi = dms_to_dd(p[0], p[1], p[2])
                lam = dms_to_dd(p[3], p[4], p[5])
            else:
                phi, lam = p[0], p[1]
            if not (-90 <= phi <= 90 and -180 <= lam <= 180):
                raise ValueError(f"Współrzędne poza zakresem: {phi:.4f}, {lam:.4f}")
            if multi_point_warning:
                QMessageBox.warning(self, "Informacja", "Otwarto tylko pierwszy punkt z pliku.")
            url = f"https://www.google.com/maps?q={phi},{lam}"
            webbrowser.open(url)
            self.show_info(f"Otworzono: {phi:.6f}°, {lam:.6f}°")
        except Exception as e:
            self.status_area.setText(f"BŁĄD\n{str(e)}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())