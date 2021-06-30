# Intproc

Одним из методов исследования высшей нервной деятельности является флюоресцентная микроскопия. Это связано с тем, что самые распространённые глиальные клетки – астроциты хоть и являются электрически невозбудимыми, но также участвуют в процессах, связанных с передачей сигналов, генерируя кальциевые сигналы. А кальциевые сигналы можно легко визуализировать при помощи флуоресцентных маркеров.

Объём данных, полученных при флюоресцентной микроскопии очень велик как по количеству изучаемых клеток, так и по времени наблюдения, поэтому актуальна задача алгоритмической обработки изображений для получения некоторых численных метрик. Для этого могут быть применены алгоритмы машинного обучения, осуществляющие классификацию объектов и событий.

Данная библиотека содержит функции, обеспечивающие обработку временного ряда интенсивности

## Установка / Начало работы

Можно производить работу в jupiter notebook

Скопируйте на свой диск модуль intproc.py и пример test.ipynb

В файле test.ipynb в строке sys.path.append(&#39;C:/Users/1/!myalgorithm&#39;) указать свой путь к модулю intproc.py

## Разработка

В файле test.ipynb - примеры применения функций модуля intproc.py и пример нахождения границ событий с использованием скользящих средних.
 В файле find borders.ipynb - попытка применения метода решающих деревьев для нахождения границ событий.

## Развертывание / публикация

intproc.py - модуль обработки временного ряда интенсивности.

## Функции

intproc.py - модуль обработки временного ряда интенсивности.

list\_files(directory)
 Список csv файлов в каталоге
 directory - путь к каталогу
 Возвращает 2 списка:
 lp - полный путь к файлу
 ln - только имена файлов

csv\_to\_lists(fn)
 Функция чтения csv файла интенсивности
 загружает как исходный файл интенсивности (с одним столбцом),
 так и размеченный в программе разметки (с тремя)
 Входные параметры:
 fn - путь к файлу
 Возвращает 3 списка
 Intensity - временной ряд интенсивности
 Borders - разметка границ
 Type - тип события

find\_max(MasDat, MA)
 Функция поиска локальных максимумов
 Входные параметры:
 MasDat - данные временного ряда, в котором требуется найти максимумы
 MA - средняя линия
 Возвращает
 LocalMaxIndex - индексы локальных максимумов
 LocalMax - значения локальных максимумов

find\_min(MasDat, MA)
 Функция поиска локальных минимумов
 Входные параметры:
 MasDat - данные временного ряда, в котором требуется найти максимумы
 MA - средняя линия
 Возвращает
 LocalMinIndex - индексы локальных минимумов
 LocalMin - значения локальных минимумов

find\_baseline(MasDat, ListMin, ListMinIndex)
 Функция поиска базовой линии
 Входные параметры:
 MasDat - данные временного ряда, в котором требуется найти максимумы
 ListMin - значения локальных минимумов
 ListMinIndex - индексы локальных минимумов
 Возвращает
 Baseline - базовая линия

find\_MA(MasDat, koefMA) Функция поиска средней линии
 Входные параметры:
 MasDat - данные временного ряда, для которого находится средняя линия
 koefMA - коэффициент в диапазоне 0..1, чем ближе к 1, тем сильнее фильтрация
 по умолчанию 0.9
 Возвращает
 MA - средняя линия

find\_MAb(MasDat, koefMA = 0.9)
 Функция поиска средней линии (с конца)
 Входные параметры:
 MasDat - данные временного ряда, для которого находится средняя линия
 koefMA - коэффициент в диапазоне 0..1, чем ближе к 1, тем сильнее фильтрация
 по умолчанию 0.9
 Возвращает
 MA - средняя линия

calc\_perc(MasDat, Baseline, ListMax)
 Функция расчёта параметров локальных максимумов
 Входные параметры:
 MasDat - данные временного ряда
 Baseline - базовая линия
 ListMax - индексы локальных максимумов
 Возвращает
 listAbs - абсолютные расстояния от линии интенсивности
 listPerc - отклонение максимумов в процентах от базовой линии

calc\_borders(Intensity):
 Функция поиска границ событий
 Входные параметры:
 Intensity - данные временного ряда интенсивности
 Возвращает
 borders - список границ событий [[начало1, максимум1, конец1], ... [началоN, максимумN, конецN]]

file\_to\_events(filename)
 Функция поиска событий
 загружает как исходный файл интенсивности (с одним столбцом),
 так и размеченный в программе разметки (с тремя)
 Входные параметры:
 filename - путь к файлу
 Возвращает:
 LocalMaxIndex - индексы локальных максимумов
 listAbs - абсолютные расстояния от линии интенсивности
 listPerc - отклонение максимумов в процентах от базовой линии
 borders - список границ событий автоматически

cut\_events(Intencity, Borders, Typ, delta = 0)
 Функция вырезания событий из списков, прочитанных функцией csv\_to\_lists из файла разметки, созданного программой разметки.
 Intensity - временной ряд интенсивности
 Borders - разметка границ
 Type - тип события
 delta - сколько точек захватить слева и справа от границы события, по умолчанию 0.
 Возвращает:
 listevents - содержит фрагменты временного ряда интенсивности, соответствующие событиям
 lborbers - содержит разметку границ


## Содействие

Если вы хотите внести свой вклад, пожалуйста, создайте ответвление репозитория и используйте функциональную ветку. Запросы на извлечение приветствуются

## Ссылки

- Домашняя страница проекта: [https://github.com/TVK-dev/Intensity](https://github.com/TVK-dev/Intensity) 
- Репозиторий: [https://github.com/TVK-dev/Intensity](https://github.com/TVK-dev/Intensity)
- Связанные проекты:
  - Программа разметки: https://github.com/TVK-dev/TruEvent
  - Репозиторий «AstrocyteLaboratory»: [https://github.com/UNN-VMK-Software/astro-analysis](https://github.com/UNN-VMK-Software/astro-analysis)

## Лицензирование

Код в этом проекте находится под лицензией Attribution-NonCommercial 2.0 Generic
