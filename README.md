# Geodetic Calculator Suite

Kompleksowa aplikacja desktopowa do wykonywania obliczeń geodezyjnych, transformacji współrzędnych oraz analiz przestrzennych na elipsoidzie odniesienia.

Opis projektu

Geodetic Calculator Suite to narzędzie umożliwiające wykonywanie najważniejszych obliczeń stosowanych w geodezji wyższej i geodezji satelitarnej. Program obsługuje zarówno układy geocentryczne, geodezyjne oraz państwowe układy współrzędnych stosowane w Polsce (PL-1992 i PL-2000).

Aplikacja oferuje intuicyjny interfejs graficzny oparty o Qt (PySide6), umożliwiający wykonywanie obliczeń na pojedynczych punktach oraz przetwarzanie danych wsadowych z plików tekstowych.

Główne funkcjonalności
Konwersje współrzędnych
Transformacja współrzędnych geocentrycznych XYZ ↔ φλh (flh).
Obsługa różnych elipsoid odniesienia:
GRS80 / WGS84
Krasowskiego
własna elipsoida definiowana przez użytkownika.
Zadanie geodezyjne wprost
Wyznaczanie współrzędnych punktu końcowego na podstawie:
współrzędnych punktu początkowego,
azymutu,
długości linii geodezyjnej.
Implementacja algorytmu Kivioja dla geodezyjnej linii na elipsoidzie.
Zadanie geodezyjne wstecz
Obliczanie:
odległości geodezyjnej,
azymutu w przód,
azymutu wstecz.
Implementacja algorytmu Vincenty'ego dla elipsoidy obrotowej.
Układy współrzędnych PL-1992 i PL-2000
Transformacje pomiędzy:
XYZ,
φλh,
PL-1992,
PL-2000.
Automatyczne wykrywanie strefy układu PL-2000.
Transformacje oparte na odwzorowaniu Gaussa-Krügera.
Redukcje geodezyjne
Skala odwzorowania.
Zbieżność południków.
Redukcja długości.
Redukcja azymutów.
Obliczenia dla układów PL-1992 i PL-2000.
Transformacje przestrzenne
Transformacje współrzędnych 3D.
Obsługa parametrów translacji, rotacji i skali.
Implementacja transformacji podobieństwa w przestrzeni.
Obliczanie pól powierzchni
Wyznaczanie pola wielokątów bezpośrednio na elipsoidzie.
Wykorzystanie biblioteki GeographicLib.
Uwzględnienie rzeczywistego kształtu Ziemi zamiast aproksymacji płaskiej.
Integracja z Google Maps
Wizualizacja punktów geograficznych w Google Maps z poziomu aplikacji.
Technologie
Python 3
PySide6 (Qt)
NumPy
GeographicLib
OOP (Object-Oriented Programming)
Zastosowania
projekty studenckie z geodezji i kartografii,
ćwiczenia laboratoryjne,
nauka transformacji współrzędnych,
analiza danych GNSS,
obliczenia geodezyjne na elipsoidzie,
weryfikacja wyników pomiarów terenowych.
Najważniejsze algorytmy
Vincenty (zadanie odwrotne)
Kivioj (zadanie wprost)
Gauss-Krüger
Transformacje PL-1992 / PL-2000
Transformacje geocentryczne XYZ ↔ φλh
Obliczanie pól geodezyjnych na elipsoidzie
Transformacje przestrzenne 3D
