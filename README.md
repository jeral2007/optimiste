# Дробные заряды через добавку к псевдопотенциалам

## python модули

+ parse\_coord.py  -- содержит функции для работы с файлами-шаблонами для координат атомов и используемых базисов

использование:
        parse_coord.py coord_pat basis_pat

вывод:

добавка к полной энергии, добавка к градиентам атомов. 

!!!Примечания:
    - *Файлы coord и basis* перезаписываются в результате вызова этого скрипта
    - требуется также модификация файла *control* для поддержки разных типов атомов (см. примеры)

+ read\_control.py -- работа с control файлом, вывод исправленной полной энергии и градиентов атомов на основе данных в файле шаблоне.

использование:
        read.control.py coord_pat control [excl_list]

*excl\_list* список имен атомов, которые нужно исключить из расчета градиентов.


## примеры использования

В папке *tests* находятся примеры файлы-шаблоны с координатами атомов для кластера из 53 атомов  и двухатомной молекулы (подкаталог *diatomic*)

Для двухатомной молекулы в файле *сontrol\_parted* приведены результаты расчета программой с патчем Ю.Д.

