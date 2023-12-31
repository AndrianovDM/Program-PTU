import pandas as pd

def breakdown_table(array, method = True):

    if method == 'param_1':
        
        data_1 = {'Температура на входе в цилиндр':'℃',
                  'Давление на входе в цилиндр':'МПа',
                  'Энтальпия на входе в цилиндр':'кДж/кг',
                  'Энтропия на входе в цилиндр':'кДж/(кгК)',
                  'Удельный объем на входе в цилиндр':'м3/кг',
                  'Степень влажности на входе в цилиндр':'-',
                  'Теоретическая температура на выходе из цилиндр':'℃',
                  'Теоретическое давление на выходе из цилиндр':'МПа',
                  'Теоретическая энтальпия на выходе из цилиндр':'кДж/кг',
                  'Теоретическая энтропия на выходе из цилиндр':'кДж/(кгК)',
                  'Теоретический удельный объем на выходе из цилиндр':'м3/кг',
                  'Теоретическая степень влажности на выходе из цилиндр':'-',
                  'Дейтсивтельная температура на выходе из цилиндр':'℃',
                  'Дейтсивтельное давление на выходе из цилиндр':'МПа',
                  'Дейтсивтельная энтальпия на выходе из цилиндр':'кДж/кг',
                  'Дейтсивтельная энтропия на выходе из цилиндр':'кДж/(кгК)',
                  'Дейтсивтельный удельный объем на выходе из цилиндр':'м3/кг',
                  'Дейтсивтельная степень влажности на выходе из цилиндр':'-',
                  'Расход пара в первую ступень': 'кг/c',
                  'Коэффициент расхода для СР': '-',
                  'Располагаемый теплоперепад на первую ступень': 'кДж/кг',
                  'Отношение скоростей U/Сф первой ступени': '-',
                  'Степень реактивности первой ступени': '-',
                  'Коэффициент потери скорости СР в первой ступени': '-',
                  'Абсолютный угол выхода потока из СР первой ступени': 'град',
                  'Окружная скорость РР первой ступени': 'м/c',
                  'Располагаемый теплоперепад на последнюю ступень': 'кДж/кг',
                  'Отношение скоростей U/Сф последней ступени': '-',
                  'Степень реактивности последней ступени': '-',
                  'Коэффициент потери скорости СР в последней ступени': '-',
                  'Абсолютный угол выхода потока из СР последней ступени': 'град',
                  'Окружная скорость РР последней ступени': 'м/c',
                  'Высота CР первой ступени': 'м',
                  'Корневой диаметр CР первой ступени': 'м',
                  'Средний диаметр CР первой ступени': 'м',
                  'Переферийный диаметр CР первой ступени': 'м',
                  'Высота РР первой ступени': 'м',
                  'Корневой диаметр РР первой ступени': 'м',
                  'Средний диаметр РР первой ступени': 'м',
                  'Переферийный диаметр РР первой ступени': 'м',
                  'Высота РР последней ступени': 'м',
                  'Корневой диаметр РР последней ступени': 'м',
                  'Средний диаметр РР последней ступени': 'м',
                  'Переферийный диаметр РР последней ступени': 'м',
                  'Сумма распологаемых теплоперепадов всех ступеней':'кДж/кг',
                  'Средний распологаемый теплоперепад на ступень':'кДж/кг',
                  'Коэффициент возврата теплоты': '-',
                  'Количество ступеней': 'шт',
                  'Величина невязки на ступень': 'кДж/кг',
                  'Теоретический теплоперепад на цилиндр': 'кДж/кг',
                  'Действительный теплоперепад на цилиндр': 'кДж/кг',
                  'Внутренний КПД цилиндра': '-',
                  'Частота вращения': 'Гц'}
        pd.set_option("display.precision", 3)
        df_merged = pd.DataFrame(list(data_1.items()), columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right') 
    
    if method == 'param_2':

        data_1 = {'Номер ступени': 'number_of_steps_i',
                'Высота РР': 'height_2_i', 
                'Корневой диаметр РР': 'D_k_2_i',
                'Средний диаметр РР': 'D_sr_2_i',
                'Переферийный диаметр РР': 'D_p_2_i',
                'Обратная веерность': 'tetta_i',
                'Степень реактивности в корне':'rho_k_i',           
                'Степень реактивности на Dsr':'rho_i',
                'Отношение скоростей U/Сф': 'X_i',
                'Коэффициент потери скорости в СР': 'fi_i',
                'Абсолютный угол выхода потока из СР':'alfa_1_i',
                'Располагаемый теплоперепад':'H_0_i',
                'Располагаемый теплоперепад с учетом невязок':'H_0_i_new',}

        df_1 = pd.DataFrame(list(data_1.items()),columns=['Параметры','Обозначение']) 
        df_2_ = pd.DataFrame(['-', 'м','м','м','м', '-', '-', '-', '-', '-', 'град', 'кДж/кг', 'кДж/кг'], columns = ['Разм'])
        df_2 = pd.DataFrame((array))
        df12_merged = df_1.join(df_2_, rsuffix ='_right')
        df_merged = df12_merged.join(df_2, rsuffix ='_right')
        pd.set_option("display.precision", 3)
    return df_merged
