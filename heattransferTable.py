import pandas as pd

def heattransfer_table(array):
    data_1 = {'Теоретическая температура на входе в СР i-ой ступени':'t_0_t_i_',
              'Теоретическое давление на входе в СР i-ой ступени':'P_0_t_i_',
              'Теоретическая энтальпия на входе в СР i-ой ступени':'h_0_t_i_',
              'Теоретическая энтропия на входе в СР i-ой ступени':'S_0_t_i_',                 
              'Теоретический объем на входе в СР i-ой ступени':'V_0_t_i_',   
              'Теоретическая степень сухости на входе в СР i-ой ступени':'x_0_t_i_',
              'Теоретическая температура на выходе из РР i-ой ступени':'t_2_t_i_',
              'Теоретическое давление на выходе из РР i-ой ступени':'P_2_t_i_',
              'Теоретическая энтальпия на выходе из РР i-ой ступени':'h_2_t_i_',
              'Теоретическая энтропия на выходе из РР i-ой ступени':'S_2_t_i_',                 
              'Теоретический объем на выходе из РР i-ой ступени':'V_2_t_i_',   
              'Теоретическая степень сухости на выходе из РР i-ой ступени':'x_2_t_i_', 
              'Действительная температура на входе в СР i-ой ступени':'t_0_i',
              'Действительная давление на входе в СР i-ой ступени':'P_0_i',
              'Действительная энтальпия на входе в СР i-ой ступени':'h_0_i',
              'Действительная энтропия на входе в СР i-ой ступени':'S_0_i',                 
              'Действительная объем на входе в СР i-ой ступени':'V_0_i',   
              'Действительная степень сухости на входе в СР i-ой ступени':'x_0_i',
              'Действительная температура на выходе из РР i-ой ступени':'t_2_i',
              'Действительная давление на выходе из РР i-ой ступени':'P_2_i',
              'Действительная энтальпия на выходе из РР i-ой ступени':'h_2_i',
              'Действительная энтропия на выходе из РР i-ой ступени':'S_2_i',                 
              'Действительная объем на выходе из РР i-ой ступени':'V_2_i',   
              'Действительная степень сухости на выходе из РР i-ой ступени':'x_2_i',             
              'Температура (заторможенного потока) на входе в СР i-ой ступени':'t_0_i_',
              'Давление (заторможенного потока) на входе в СР i-ой ступени':'P_0_i_',
              'Энтальпия (заторможенного потока) на входе в СР i-ой ступени':'h_0_i_',
              'Энтропия (заторможенного потока) на входе в СР i-ой ступени':'S_0_i_',                 
              'Объем (заторможенного потока) на входе в СР i-ой ступени':'V_0_i_',   
              'Степень сухости (заторможенного потока) на входе в СР i-ой ступени':'x_0_i_',                          
              'Температура по изоэнтропе на выходе из РР i-ой ступени':'t_2_i_',
              'Давление по изоэнтропе на выходе из РР i-ой ступени':'P_2_i_',
              'Энтальпия по изоэнтропе на выходе из РР i-ой ступени':'h_2_i_',
              'Энтропия по изоэнтропе на выходе из РР i-ой ступени':'S_2_i_',                 
              'Объем по изоэнтропе на выходе из РР i-ой ступени':'V_2_i_',   
              'Степень сухости по изоэнтропе на выходе из РР i-ой ступени':'x_2_i_',
              'Скорость на входе i-ой ступени':'C_0_i',
              'Располагаемый теплоперепад от параметров полного томожения i-ой ступени':'H_0_i_',
              'Располагаемый теплоперепад от статических параметров i-ой ступени':'H_0_i',
              }

    df_1 = pd.DataFrame(list(data_1.items()),columns=['Параметры','Обозначение']) 
    df_2_ = pd.DataFrame(['℃', 'МПа', 'кДж/кг', 'кДж/(кгК)', 'м3/кг', '-', 
                          '℃', 'МПа', 'кДж/кг', 'кДж/(кгК)', 'м3/кг', '-', 
                          '℃', 'МПа', 'кДж/кг', 'кДж/(кгК)', 'м3/кг', '-',
                          '℃', 'МПа', 'кДж/кг', 'кДж/(кгК)', 'м3/кг', '-', 
                          '℃', 'МПа', 'кДж/кг', 'кДж/(кгК)', 'м3/кг', '-', 
                          '℃', 'МПа', 'кДж/кг', 'кДж/(кгК)', 'м3/кг', '-',
                          'м/с', 'кДж/кг', 'кДж/кг'], columns = ['Разм'])
    
    df_2 = pd.DataFrame((array))
    df12_merged = df_1.join(df_2_, rsuffix='_right')
    df_merged = df12_merged.join(df_2, rsuffix='_right')
    pd.set_option("display.precision", 3)
    return df_merged