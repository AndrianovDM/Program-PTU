import pandas as pd
def stage_table(array, method = True):
    if method == 'parameters':
        data_1 = {'Температура на входе в ступень (параметры полного торможения)':'℃',
                'Температура на входе в ступень (статические параметры)':'℃',
                'Температура на выходе из СР (теоретическая)':'℃',
                'Температура на выходе из СР (действительная)':'℃',
                'Температура на выходе из СР (параметры полного торможения)':'℃',
                'Температура на выходе из РР (теоретическая)':'℃',
                'Температура на выходе из РР (статические параметры)':'℃',
                'Температура на выходе из РР (действительная)':'℃',
                'Температура на выходе из РР (параметры полного торможения)':'℃',

                'Давление на входе в ступень (параметры полного торможения)':'МПа',
                'Давление на входе в ступень (статические параметры)':'МПа',
                'Давление на выходе из СР (теоретическая)':'МПа',
                'Давление на выходе из СР (действительная)':'МПа',
                'Давление на выходе из СР (параметры полного торможения)':'МПа',
                'Давление на выходе из РР (теоретическая)':'МПа',
                'Давление на выходе из РР (статические параметры)':'МПа',
                'Давление на выходе из РР (действительная)':'МПа',
                'Давление на выходе из РР (параметры полного торможения)':'МПа',
                
                'Энтальпия на входе в ступень (параметры полного торможения)':'кДж/кг',
                'Энтальпия на входе в ступень (статические параметры)':'кДж/кг',
                'Энтальпия на выходе из СР (теоретическая)':'кДж/кг',
                'Энтальпия на выходе из СР (действительная)':'кДж/кг',
                'Энтальпия на выходе из СР (параметры полного торможения)':'кДж/кг',
                'Энтальпия на выходе из РР (теоретическая)':'кДж/кг',
                'Энтальпия на выходе из РР (статические параметры)':'кДж/кг',
                'Энтальпия на выходе из РР (действительная)':'кДж/кг',
                'Энтальпия на выходе из РР (параметры полного торможения)':'кДж/кг',
                
                'Объем на входе в ступень (параметры полного торможения)':'м3/кг',
                'Объем на входе в ступень (статические параметры)':'м3/кг',
                'Объем на выходе из СР (теоретическая)':'м3/кг',
                'Объем на выходе из СР (действительная)':'м3/кг',
                'Объем на выходе из СР (параметры полного торможения)':'м3/кг',
                'Объем на выходе из РР (теоретическая)':'м3/кг',
                'Объем на выходе из РР (статические параметры)':'м3/кг',
                'Объем на выходе из РР (действительная)':'м3/кг',
                'Объем на выходе из РР (параметры полного торможения)':'м3/кг',
                
                'Степень сухости на входе в ступень (параметры полного торможения)':'-',
                'Степень сухости на входе в ступень (статические параметры)':'-',
                'Степень сухости на выходе из СР (теоретическая)':'-',
                'Степень сухости на выходе из СР (действительная)':'-',
                'Степень сухости на выходе из СР (параметры полного торможения)':'-',
                'Степень сухости на выходе из РР (теоретическая)':'-',
                'Степень сухости на выходе из РР (статические параметры)':'-',
                'Степень сухости на выходе из РР (действительная)':'-',
                'Степень сухости на выходе из РР (параметры полного торможения)':'-',

                'Изоэнтропийный теплоперепад ступени (от параметров полного торможения)':'кДж/кг',
                'Изоэнтропийный теплоперепад СР':'кДж/кг',
                'Изоэнтропийный теплоперепад РР':'кДж/кг',
                'Использованный теплоперепад ступени':'кДж/кг',
                
                'Фиктивная скорость':'м/c',
                'Скорость на входе в ступень':'м/c',
                'Абсолютная скорость на выходе из СР':'м/c',
                'Абсолютная скорость на выходе из РР':'м/c',
                'Относительная скорость на выходе из СР':'м/c',
                'Относительная скорость на выходе из РР':'м/c',
                'Окружная скорость на выходе из СР':'м/c',
                'Окружная скорость на выходе из РР':'м/c',
                'Угол абсолютной скорости на выходе из СР':'град',
                'Угол абсолютной скорости на выходе из РР':'град',
                'Угол относительной скорости на выходе из СР':'град',
                'Угол относительной скорости на выходе из РР':'град',

                'Мах по асболютной скорости из СР':'-',
                'Мах по асболютной скорости из РР':'-',
                'Мах по относительной скорости из СР':'-',
                'Мах по относительной скорости из РР':'-',

                'Коэффициент потери скорости в СР':'-',
                'Коэффициент потери скорости в РР':'-',
                'Коэффициент суммарных потерь энергии в СР':'-',
                'Коэффициент суммарных потерь энергии в РР':'-',
                'Коэффициент профильныех потерь энергии в СР':'-',
                'Коэффициент профильныех потерь энергии в РР':'-',
                'Коэффициент вторичных потерь энергии в СР':'-',
                'Коэффициент вторичных потерь энергии в РР':'-',
                'Коэффициент потерь энергии связанные с утечеками в СР':'-',
                'Коэффициент потерь энергии связанные с утечеками в РР':'-',
                'Коэффициент кромочных потерь энергии в СР':'-',
                'Коэффициент кромочных потерь энергии в РР':'-',
                'Коэффициент потерь энергии связанные с утечеками через переферию в СР':'-',
                'Коэффициент потерь энергии связанные с утечеками через переферию в РР':'-',
                'Высота СР':'м',
                'Высота РР':'м',
                'Средний диаметр СР':'м',
                'Средний диаметр РР':'м', 
                'Выходная площадь СР':'м2',
                'Выходная площадь РР':'м2',
                'Коэффициент расхода СР':'-',
                'Коэффициент расхода РР':'-',
                'Потери в СР':'кДж/кг',
                'Потери в РР':'кДж/кг',
                'Потери от трения диска':'кДж/кг',
                'Потери от утечек':'кДж/кг',
                'Потери от влажности':'кДж/кг',
                'Энергия выходной скорости':'кДж/кг',
                'Располагаемая энергия ступени':'кДж/кг',
                'Окружная работа ступени':'кДж/кг',
                'Степень реактивности':'-',
                'Отношение скоростей U/Сф':'-',
                'Оптимальное отношение скоростей U/Сф':'-',
                'Относительный лопаточный КПД ступени':'-',
                'Внутренний относительный КПД ступени':'-',
                'Расход рабочего тела на ступень':'кг/c',
                'Внутренняя мощность ступени':'кВт'}
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  

    if method == 'profile sopl':
        data_1 = {'Ширина профиля СА':'мм', 
                  'Хорда профиля СА': 'мм', 
                  'Шаг решетки СА':'мм',
                  'Радиус входной кромки СА':'мм', 
                  'Радиус выходной кромки СА': 'мм',
                  'Горло решетки на входе СА':'мм',
                  'Горло решетки на выходе СА':'мм', 
                  'Координата центра СА':'мм', 
                  'Толщина профиля СА':'мм',
                  'Корневой диаметр СА':'мм', 
                  'Высота лопатки СА':'мм', 
                  'Угол установки СА':'град',
                  'Лопаточный угол на входе СА':'град',
                  'Лопаточный угол на выходе СА':'град',
                  'Угол заострения на вх. СА':'град', 
                  'Угол заострения на вых. СА':'град',
                  'Угол отгиба СА':'град',
                  'Количество лопаток СА':'шт'}
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  

    if method == 'profile rab':
        data_1 = {'Ширина профиля РК':'мм', 
                  'Хорда профиля РК': 'мм', 
                  'Шаг решетки РК':'мм',
                  'Радиус входной кромки РК':'мм', 
                  'Радиус выходной кромки РК': 'мм',
                  'Горло решетки на входе РК':'мм',
                  'Горло решетки на выходе РК':'мм', 
                  'Координата центра РК':'мм',
                  'Толщина профиля РК':'мм',
                  'Корневой диаметр РК':'мм', 
                  'Высота лопатки РК':'мм', 
                  'Угол установки РК':'град',
                  'Лопаточный угол на входе РК':'град',
                  'Лопаточный угол на выходе РК':'град',
                  'Угол заострения на вх. РК':'град', 
                  'Угол заострения на вых. РК':'град',
                  'Угол отгиба РК':'град',
                  'Количество лопаток РК':'шт'}
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  
    
    pd.set_option("display.precision", 3
    
    )
    return df_merged   

def stageTable(array):
    data_1 = {'Температура на входе в ступень (параметры полного торможения)':'t_0_',
            'Температура на входе в ступень (статические параметры)':'t_0',
            'Температура на выходе из СР (теоретическая)':'t_1t',
            'Температура на выходе из СР (действительная)':'t_1',
            'Температура на выходе из СР (параметры полного торможения)':'t_1w',
            'Температура на выходе из РР (теоретическая)':'t_2t_',
            'Температура на выходе из РР (статические параметры)':'t_2t',
            'Температура на выходе из РР (действительная)':'t_2',
            'Температура на выходе из РР (параметры полного торможения)':'t_2w',

            'Давление на входе в ступень (параметры полного торможения)':'P_0_',
            'Давление на входе в ступень (статические параметры)':'P_0',
            'Давление на выходе из СР (теоретическая)':'P_1t',
            'Давление на выходе из СР (действительная)':'P_1',
            'Давление на выходе из СР (параметры полного торможения)':'P_1w',
            'Давление на выходе из РР (теоретическая)':'P_2t_',
            'Давление на выходе из РР (статические параметры)':'P_2t',
            'Давление на выходе из РР (действительная)':'P_2',
            'Давление на выходе из РР (параметры полного торможения)':'P_2w',
                
            'Энтальпия на входе в ступень (параметры полного торможения)':'h_0_',
            'Энтальпия на входе в ступень (статические параметры)':'h_0',
            'Энтальпия на выходе из СР (теоретическая)':'h_1t',
            'Энтальпия на выходе из СР (действительная)':'h_1',
            'Энтальпия на выходе из СР (параметры полного торможения)':'h_1w',
            'Энтальпия на выходе из РР (теоретическая)':'h_2t_',
            'Энтальпия на выходе из РР (статические параметры)':'h_2t',
            'Энтальпия на выходе из РР (действительная)':'h_2',
            'Энтальпия на выходе из РР (параметры полного торможения)':'h_2w',
                
            'Объем на входе в ступень (параметры полного торможения)':'V_0_',
            'Объем на входе в ступень (статические параметры)':'V_0',
            'Объем на выходе из СР (теоретическая)':'V_1t',
            'Объем на выходе из СР (действительная)':'V_1',
            'Объем на выходе из СР (параметры полного торможения)':'V_1w',
            'Объем на выходе из РР (теоретическая)':'V_2t_',
            'Объем на выходе из РР (статические параметры)':'V_2t',
            'Объем на выходе из РР (действительная)':'V_2',
            'Объем на выходе из РР (параметры полного торможения)':'V_2w',
                
            'Степень сухости на входе в ступень (параметры полного торможения)':'x_0_',
            'Степень сухости на входе в ступень (статические параметры)':'x_0',
            'Степень сухости на выходе из СР (теоретическая)':'x_1t',
            'Степень сухости на выходе из СР (действительная)':'x_1',
            'Степень сухости на выходе из СР (параметры полного торможения)':'x_1w',
            'Степень сухости на выходе из РР (теоретическая)':'x_2t_',
            'Степень сухости на выходе из РР (статические параметры)':'x_2t',
            'Степень сухости на выходе из РР (действительная)':'x_2',
            'Степень сухости на выходе из РР (параметры полного торможения)':'x_2w',

            'Изоэнтропийный теплоперепад ступени (от параметров полного торможения)':'H_0_',
            'Изоэнтропийный теплоперепад СР':'H_0_sa_',
            'Изоэнтропийный теплоперепад РР':'H_0_rk',
            'Использованный теплоперепад ступени':'H_i',
                
            'Фиктивная скорость':'C_fict',
            'Скорость на входе в ступень':'C_0',
            'Абсолютная скорость на выходе из СР':'C1',
            'Абсолютная скорость на выходе из РР':'C2',
            'Относительная скорость на выходе из СР':'W1',
            'Относительная скорость на выходе из РР':'W2',
            'Окружная скорость на выходе из СР':'U1',
            'Окружная скорость на выходе из РР':'U2',
            'Угол абсолютной скорости на выходе из СР':'alpha1',
            'Угол абсолютной скорости на выходе из РР':'alpha2',
            'Угол относительной скорости на выходе из СР':'betta1',
            'Угол относительной скорости на выходе из РР':'betta2',

            'Мах по асболютной скорости из СР':'M1c',
            'Мах по асболютной скорости из РР':'M2c',
            'Мах по относительной скорости из СР':'M1w',
            'Мах по относительной скорости из РР':'M2w',

            'Коэффициент потери скорости в СР':'fi',
            'Коэффициент потери скорости в РР':'psi',
            'Коэффициент суммарных потерь энергии в СР':'Y_s_sopl',
            'Коэффициент суммарных потерь энергии в РР':'Y_s_rab',
            'Коэффициент профильныех потерь энергии в СР':'Y_p_sopl',
            'Коэффициент профильныех потерь энергии в РР':'Y_p_rab',
            'Коэффициент вторичных потерь энергии в СР':'Y_sec_sopl',
            'Коэффициент вторичных потерь энергии в РР':'Y_sec_rab',
            'Коэффициент потерь энергии связанные с утечеками в СР':'Y_tl_sopl',
            'Коэффициент потерь энергии связанные с утечеками в РР':'Y_tl_rab',
            'Коэффициент кромочных потерь энергии в СР':'Y_te_sopl',
            'Коэффициент кромочных потерь энергии в РР':'Y_te_rab',
            'Коэффициент потерь энергии связанные с утечеками через переферию в СР':'Y_cl_sopl',
            'Коэффициент потерь энергии связанные с утечеками через переферию в РР':'Y_cl_rab',
            'Высота СР':'h_sopl',
            'Высота РР':'h_rab',
            'Средний диаметр СР':'Dsr_sopl',
            'Средний диаметр РР':'Dsr_rab', 
            'Выходная площадь СР':'F_1_out',
            'Выходная площадь РР':'F_2_out',
            'Коэффициент расхода СР':'mu_1',
            'Коэффициент расхода РР':'mu_2',
            'Потери в СР':'Delta_Hs',
            'Потери в РР':'Delta_Hr',
            'Потери от трения диска':'delta_Htr',
            'Потери от утечек':'delta_Hlake',
            'Потери от влажности':'delta_Hvet',
            'Энергия выходной скорости':'Delta_Hvs',
            'Располагаемая энергия ступени':'E_0',
            'Окружная работа ступени':'L_u',
            'Степень реактивности':'rho',
            'Отношение скоростей U/Сф':'X',
            'Оптимальное отношение скоростей U/Сф':'X_opt',
            'Относительный лопаточный КПД ступени':'etta_ol',
            'Внутренний относительный КПД ступени':'etta_oi',
            'Расход рабочего тела на ступень':'G_0',
            'Внутренняя мощность ступени':'N_i'}

    df_1 = pd.DataFrame(list(data_1.items()),columns=['Параметры','Обозначение']) 
    df_2_ = pd.DataFrame(['℃', '℃', '℃', '℃', '℃', '℃', '℃', '℃', '℃', 
                          'МПа', 'МПа', 'МПа', 'МПа', 'МПа', 'МПа', 'МПа', 'МПа', 'МПа',
                          'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг',
                          'м3/кг', 'м3/кг', 'м3/кг', 'м3/кг', 'м3/кг', 'м3/кг', 'м3/кг', 'м3/кг', 'м3/кг',
                          '-', '-', '-', '-', '-', '-', '-', '-', '-', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг',
                          'м/c', 'м/c', 'м/c', 'м/c', 'м/c', 'м/c', 'м/c', 'м/c', 'град', 'град', 'град', 'град',
                          '-', '-', '-', '-','-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
                          'м', 'м', 'м', 'м', 'м2', 'м2', '-', '-', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг',
                          '-', '-', '-', '-', '-', 'кг/c', 'кВт'], columns = ['Разм'])
    df_2 = pd.DataFrame((array))
    df12_merged = df_1.join(df_2_, rsuffix='_right')
    df_merged = df12_merged.join(df_2, rsuffix='_right')
    pd.set_option("display.precision", 3)
    return df_merged