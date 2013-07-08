## Использование программы сегментации изображения КТ лёгкого

На вход ожидается файл с изображением КТ обследования в формате Analyze 7.5

Результатом работы сегментатора будут два файла в том же формате. Один представляет собой маску для сегментации с суффиксом _mask, во втором файле _segm содержится сегментированное изображение лёгкого. 

> Запуск сегментации:
>   
> 	ct-segmentation.exe W0001.hdr
> 
> Результат сегментации: 
> 
> 	W0001_mask.hdr   
> 	W0001_mask.img  
> 	W0001_segm.hdr  
> 	W0001_segm.img

### На примере одного изображения

Слой №138 из серии W0001.
![Оригинальное изображение](./img/W0001-138-original.jpg "Оригинальное изображение")

Область лёгкого отмечена нулями, фон - единицами.
![Маска изображения](./img/W0001-138-mask.jpg "Маска изображения")

Пикселам фона присваивается наименьшее значение из найденных на оригинальном изображении.
![Результат сегментации](./img/W0001-138-segm.jpg "Результат сегментации")

## Алгоритм сегментации
Вся серия сегментируется послойно. Для каждого слоя:

- находится пороговое значение, отделяющее область, похожую на лёгкое, от области, не похожей на лёгкое;
- слой сегментируется по найденному пороговому значению;
- сегментированные участки удаляются либо сливаются с другими;
- финальная маска применяется к исходному слою, и получается один слой сегментированного лёгкого.

## Изображения и сторонние библиотеки

Источник тестовых данных: [ELCAP Public Lung Image Database](http://www.via.cornell.edu/lungdb.html)
Оригинальные DICOM серии конвертировались в формат Analyze 7.5 с помощью программы [MRIConvert v 2.0 rev. 235](http://lcni.uoregon.edu/~jolinda/MRIConvert/).

Библиотека обработки изображений: [CImg v 1.5.6](http://cimg.sourceforge.net)











































