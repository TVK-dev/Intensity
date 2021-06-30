
В файле test.ipynb - примеры применения функций модуля intproc.py и пример нахождения границ событий с использованием скользящих средних.  
В файле find borders.ipynb - попытка применения метода решающих деревьев для нахождения границ событий.  

intproc.py - модуль обработки временного ряда интенсивности.  
Функции модуля обработки временного ряда интенсивности 
Intensity processing (intproc.py):  


list_files(directory)  
Список csv файлов в каталоге  
directory - путь к каталогу  
Возвращает 2 списка:  
lp - полный путь к файлу  
ln - только имена файлов  
  
  
csv_to_lists(fn)  
Функция чтения csv файла интенсивности  
загружает как исходный файл интенсивности (с одним столбцом),  
так и размеченный в программе разметки (с тремя)  
Входные параметры:  
fn - путь к файлу  
Возвращает 3 списка  
Intensity - временной ряд интенсивности  
Borders - разметка границ  
Type - тип события  



find_max(MasDat, MA)  
Функция поиска локальных максимумов  
Входные параметры:  
MasDat - данные временного ряда, в котором требуется найти максимумы  
MA - средняя линия  
Возвращает   
LocalMaxIndex - индексы локальных максимумов  
LocalMax - значения локальных максимумов  
         

find_min(MasDat, MA)  
Функция поиска локальных минимумов  
Входные параметры:  
MasDat - данные временного ряда, в котором требуется найти максимумы  
MA - средняя линия  
Возвращает  
LocalMinIndex - индексы локальных минимумов  
LocalMin - значения локальных минимумов  
    
    
find_baseline(MasDat, ListMin, ListMinIndex)  
Функция поиска базовой линии   
Входные параметры:  
MasDat - данные временного ряда, в котором требуется найти максимумы  
ListMin - значения локальных минимумов  
ListMinIndex - индексы локальных минимумов  
Возвращает  
Baseline - базовая линия  


find_MA(MasDat, koefMA) 
Функция поиска средней линии   
Входные параметры:  
MasDat - данные временного ряда, для которого находится средняя линия  
koefMA - коэффициент в диапазоне 0..1, чем ближе к 1, тем сильнее фильтрация  
по умолчанию 0.9  
Возвращает  
MA - средняя линия    


find_MAb(MasDat, koefMA = 0.9)  
Функция поиска средней линии (с конца)  
Входные параметры:  
MasDat - данные временного ряда, для которого находится средняя линия  
koefMA - коэффициент в диапазоне 0..1, чем ближе к 1, тем сильнее фильтрация  
по умолчанию 0.9  
Возвращает  
MA - средняя линия    


calc_perc(MasDat, Baseline, ListMax)  
Функция расчёта параметров локальных максимумов  
Входные параметры:  
MasDat - данные временного ряда  
Baseline - базовая линия  
ListMax - индексы локальных максимумов  
Возвращает   
listAbs - абсолютные расстояния  от линии интенсивности  
listPerc - отклонение максимумов в процентах от базовой линии  


calc_borders(Intensity):  
Функция поиска границ событий  
Входные параметры:  
Intensity - данные временного ряда интенсивности  
Возвращает   
borders - список границ событий [[начало1, максимум1, конец1], ... [началоN, максимумN, конецN]]  


file_to_events(filename)  
Функция поиска событий  
загружает как исходный файл интенсивности (с одним столбцом),  
так и размеченный в программе разметки (с тремя)  
Входные параметры:  
filename - путь к файлу  
Возвращает:  
LocalMaxIndex - индексы локальных максимумов  
listAbs - абсолютные расстояния  от линии интенсивности  
listPerc - отклонение максимумов в процентах от базовой линии  
borders - список границ событий автоматически  


cut_events(Intencity, Borders, Typ, delta = 0)  
Функция вырезания событий из списков, прочитанных функцией csv_to_lists из файла разметки, созданного программой разметки.
Intensity - временной ряд интенсивности  
Borders - разметка границ  
Type - тип события  
delta - сколько точек захватить слева и справа от границы события, по умолчанию 0
Возвращает:  
listevents - содержит фрагменты временного ряда интенсивности, соответствующие событиям
lborbers - содержит разметку границ
    




