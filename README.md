# Geodetic Calculator Suite (ENG BELOW)

Kompleksowa aplikacja desktopowa do wykonywania obliczeń geodezyjnych, transformacji współrzędnych oraz analiz przestrzennych na elipsoidzie odniesienia.

## Opis projektu

Geodetic Calculator Suite to narzędzie umożliwiające wykonywanie najważniejszych obliczeń stosowanych w geodezji wyższej i geodezji satelitarnej. Program obsługuje układy geocentryczne, geodezyjne oraz państwowe układy współrzędnych stosowane w Polsce (PL-1992 i PL-2000).

Aplikacja oferuje intuicyjny interfejs graficzny oparty na Qt (PySide6), umożliwiający wykonywanie obliczeń dla pojedynczych punktów oraz przetwarzanie danych wsadowych z plików tekstowych.

## Główne funkcjonalności

### Konwersje współrzędnych

* Transformacja współrzędnych geocentrycznych **XYZ ↔ φλh (flh)**.
* Obsługa różnych elipsoid odniesienia:

  * GRS80 / WGS84
  * Krasowskiego
  * Własna elipsoida definiowana przez użytkownika

### Zadanie geodezyjne wprost

Wyznaczanie współrzędnych punktu końcowego na podstawie:

* współrzędnych punktu początkowego,
* azymutu,
* długości linii geodezyjnej.

Implementacja algorytmu **Kivioja** dla linii geodezyjnej na elipsoidzie.

### Zadanie geodezyjne wstecz

Obliczanie:

* odległości geodezyjnej,
* azymutu w przód,
* azymutu wstecz.

Implementacja algorytmu **Vincenty'ego** dla elipsoidy obrotowej.

### Układy współrzędnych PL-1992 i PL-2000

Transformacje pomiędzy:

* XYZ,
* φλh,
* PL-1992,
* PL-2000.

Dodatkowe funkcjonalności:

* automatyczne wykrywanie strefy układu PL-2000,
* transformacje oparte na odwzorowaniu **Gaussa–Krügera**.

### Redukcje geodezyjne

* skala odwzorowania,
* zbieżność południków,
* redukcja długości,
* redukcja azymutów,
* obliczenia dla układów PL-1992 i PL-2000.

### Transformacje przestrzenne

* transformacje współrzędnych 3D,
* obsługa parametrów translacji, rotacji i skali,
* implementacja transformacji podobieństwa w przestrzeni.

### Obliczanie pól powierzchni

* wyznaczanie pola wielokątów bezpośrednio na elipsoidzie,
* wykorzystanie biblioteki GeographicLib,
* uwzględnienie rzeczywistego kształtu Ziemi zamiast aproksymacji płaskiej.

### Integracja z Google Maps

* wizualizacja punktów geograficznych bezpośrednio w Google Maps z poziomu aplikacji.

## Technologie

* Python 3
* PySide6 (Qt)
* NumPy
* GeographicLib
* Object-Oriented Programming (OOP)

## Zastosowania

* projekty studenckie z geodezji i kartografii,
* ćwiczenia laboratoryjne,
* nauka transformacji współrzędnych,
* analiza danych GNSS,
* obliczenia geodezyjne na elipsoidzie,
* weryfikacja wyników pomiarów terenowych.

## Zaimplementowane algorytmy

* Vincenty (zadanie geodezyjne wstecz)
* Kivioj (zadanie geodezyjne wprost)
* Gauss–Krüger
* Transformacje PL-1992 / PL-2000
* Transformacje geocentryczne XYZ ↔ φλh
* Obliczanie pól geodezyjnych na elipsoidzie
* Transformacje przestrzenne 3D

---

# Geodetic Calculator Suite

A comprehensive desktop application for geodetic computations, coordinate transformations, and spatial analysis on a reference ellipsoid.

## Project Description

Geodetic Calculator Suite is a software tool designed to perform essential calculations used in higher geodesy and satellite geodesy. The application supports geocentric, geodetic, and Polish national coordinate systems (PL-1992 and PL-2000).

The program provides an intuitive graphical user interface built with Qt (PySide6), allowing users to perform calculations on individual points as well as process batch data from text files.

## Main Features

### Coordinate Transformations

* Conversion between geocentric and geodetic coordinates: **XYZ ↔ φλh (latitude, longitude, ellipsoidal height)**.
* Support for multiple reference ellipsoids:

  * GRS80 / WGS84
  * Krasowski
  * User-defined custom ellipsoid

### Direct Geodetic Problem

Determination of endpoint coordinates based on:

* starting point coordinates,
* azimuth,
* geodesic distance.

Implemented using the **Kivioja algorithm** for geodesic lines on an ellipsoid.

### Inverse Geodetic Problem

Calculation of:

* geodesic distance,
* forward azimuth,
* reverse azimuth.

Implemented using the **Vincenty algorithm** for rotational ellipsoids.

### PL-1992 and PL-2000 Coordinate Systems

Coordinate transformations between:

* XYZ,
* φλh,
* PL-1992,
* PL-2000.

Additional features:

* automatic PL-2000 zone detection,
* transformations based on the **Gauss–Krüger projection**.

### Geodetic Reductions

* projection scale factor,
* meridian convergence,
* distance reduction,
* azimuth reduction,
* calculations for PL-1992 and PL-2000 systems.

### Spatial Transformations

* 3D coordinate transformations,
* support for translation, rotation, and scale parameters,
* similarity transformation implementation.

### Surface Area Computation

* calculation of polygon areas directly on the ellipsoid,
* integration with GeographicLib,
* accurate consideration of the Earth's ellipsoidal shape instead of planar approximations.

### Google Maps Integration

* visualization of geographic points directly in Google Maps from within the application.

## Technologies

* Python 3
* PySide6 (Qt)
* NumPy
* GeographicLib
* Object-Oriented Programming (OOP)

## Applications

* geodesy and cartography student projects,
* laboratory exercises,
* learning coordinate transformation techniques,
* GNSS data analysis,
* geodetic calculations on reference ellipsoids,
* verification of field survey measurements.

## Implemented Algorithms

* Vincenty (Inverse Geodetic Problem)
* Kivioja (Direct Geodetic Problem)
* Gauss–Krüger Projection
* PL-1992 / PL-2000 Transformations
* Geocentric to Geodetic Coordinate Transformations (XYZ ↔ φλh)
* Ellipsoidal Area Computation
* 3D Spatial Transformations
