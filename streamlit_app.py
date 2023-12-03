from streamlitFunction import *
from geometry import *
from heattransfer import *
from stage import *
from section import *

st.set_option('deprecation.showPyplotGlobalUse', False)
st.set_page_config(
        page_title = "Program PTU",
        page_icon = Image.open("icon.ico"),
        layout = "wide",
        initial_sidebar_state = "expanded")

panel_global = st.sidebar.radio('Этапы расчета:', ["I. - Этап расчета ПТУ", "II. - Этап расчета ступеней ПТУ", "III. - Этап расчета по сечениям"])

if panel_global == "I. - Этап расчета ПТУ":
    panel_1 = st.sidebar.radio('Этапы расчета проточной части:', 
                               ["1. - Расчет цилиндров ПТУ",  
                                "2. - Расчет геометрии ПТУ", 
                                "3. - Расчет параметров по ступеням"])
    
    if panel_1 == "1. - Расчет цилиндров ПТУ":
        st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет цилиндров ПТУ</ins> </h1>", unsafe_allow_html=True)
        st.header('Ввод исходных данных')
       
        st.session_state.bdown = None
        st.session_state.geom = None
        st.session_state.heat = None

        select_1 = st.selectbox('Выберите цилиндр расчета ПТУ:', ('Расчет ЦВД', 'Расчет ЦСД', 'Расчет ЦНД'))
        if select_1 ==  'Расчет ЦВД':
            with st.form(key = 'my_form_1'):
                    _vP_0_, _vt_0_, _vP_2_z, = st.columns(3)
                    _vG_0, _vG_z, _vetta_oi = st.columns(3)
                    _vrho_k_1, _valfa_1_1, _vfi_1, = st.columns(3)
                    _vrho_k_z, _valfa_1_z, _vfi_z = st.columns(3)
                    _vD_sr_reg, _vdelta_D, _vdelta, = st.columns(3)
                    _vn, _= st.columns(2)
                    
                    P_0_ = count(counter = 'vP_0_', column = _vP_0_, name = "P_0, МПа", format = "%f", ke = 'count1')
                    t_0_ = count(counter = 'vt_0_', column = _vt_0_, name = "t_0, ℃", format = "%f", ke = 'count2')
                    P_2_z = count(counter = 'vP_2_z', column = _vP_2_z, name = "P_2z, МПа", format = "%f", ke = 'count3')

                    G_0 = count(counter = 'vG_0', column = _vG_0, name = "G_0, кг/с", format = "%f", ke = 'count4')
                    G_z = count(counter = 'vG_z', column = _vG_z, name = "G_z, кг/с", format = "%f", ke = 'count5')

                    rho_k_1 = count(counter = 'vrho_k_1', column = _vrho_k_1, name = "rho_k_1, -", format = "%f", ke = 'count7')
                    alfa_1_1 = count(counter = 'valfa_1_1', column = _valfa_1_1, name = "alfa_1_1, град", format = "%f", ke = 'count8')
                    fi_1 = count(counter = 'vfi_1', column = _vfi_1, name = "fi_1, -", format = "%f", ke = 'count9')

                    rho_k_z = count(counter = 'vrho_k_z', column = _vrho_k_z, name = "rho_k_z, -", format = "%f", ke = 'count10')
                    alfa_1_z = count(counter = 'valfa_1_z', column = _valfa_1_z, name = "alfa_1_z, град", format = "%f", ke = 'count11')
                    fi_z = count(counter = 'vfi_z', column = _vfi_z, name = "fi_z, -", format = "%f", ke = 'count12')

                    D_sr_reg = count(counter = 'vD_sr_reg', column = _vD_sr_reg, name = "D_sr_reg, м", format = "%f", ke = 'count13')
                    delta_D = count(counter = 'vdelta_D', column = _vdelta_D, name = "delta_D, м", format = "%f", ke = 'count14')
                    delta = count(counter = 'vdelta', column = _vdelta, name = "delta, м", format = "%f", ke = 'count15')

                    etta_oi = count(counter = 'vetta_oi', column = _vetta_oi, name = "etta_oi, -", format = "%f", ke = 'count16')
                    n = count(counter = 'vn', column = _vn, name = "n, Гц", format = "%g", ke = 'count17')
                   
                    # st.write('P_0_ = ',P_0_, '; ','t_0_ = ',t_0_, '; ','P_2_z = ',P_2_z)
                    # st.write('G_0 = ',G_0, '; ','G_z = ',G_z, '; ','mu_1 = ',mu_1)
                    # st.write('rho_k_1 = ',rho_k_1, '; ','alfa_1_1 = ',alfa_1_1, '; ','fi_1 = ',fi_1)
                    # st.write('rho_k_z = ',rho_k_z, '; ','alfa_1_z = ',alfa_1_z, '; ','fi_z = ',fi_z)
                    # st.write('D_sr_reg = ',D_sr_reg, '; ','delta_D = ',delta_D, '; ','delta = ',delta)
                    # st.write('etta_oi = ',etta_oi, '; ','n = ',n)

                    if st.form_submit_button('Расчет'):
                        bdown = breakdown_CVP(P_0_ = st.session_state.vP_0_, t_0_ = st.session_state.vt_0_, P_2_z =  st.session_state.vP_2_z, 
                        G_0 = st.session_state.vG_0, G_z = st.session_state.vG_z, etta_oi = st.session_state.vetta_oi,   
                        rho_k_1 = st.session_state.vrho_k_1, alfa_1_1 = st.session_state.valfa_1_1, fi_1 = st.session_state.vfi_1,
                        rho_k_z = st.session_state.vrho_k_z, alfa_1_z = st.session_state.valfa_1_z, fi_z = st.session_state.vfi_z, 
                        D_sr_reg = st.session_state.vD_sr_reg, delta_D = st.session_state.vdelta_D, delta = st.session_state.vdelta,
                        n = st.session_state.vn)
                        st.session_state.bdown = bdown
                        
                        st.header('Результаты расчета ЦВД')
                        st.table(breakdown_table(st.session_state.bdown[0], 'param_1'))
                        st.header('Распределение параметров по ступеням')
                        st.table(breakdown_table(st.session_state.bdown[2], 'param_2'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][4], k =  st.session_state.bdown[2][3], z = st.session_state.bdown[2][2],  method = 'Form_diametr'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][3], method = 'Diametr'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][1], method = 'Height'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][5], method = 'Tetta'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][7], method = 'Reactivity'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][8], method = 'X'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][11], method = 'Heat_difference'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][12], method = 'Heat_difference_new'))
                    
                    if st.form_submit_button('Сохранить в Excel'):
                        bdown = breakdown_CVP(P_0_ = st.session_state.vP_0_, t_0_ = st.session_state.vt_0_, P_2_z =  st.session_state.vP_2_z, 
                        G_0 = st.session_state.vG_0, G_z = st.session_state.vG_z, etta_oi = st.session_state.vetta_oi,   
                        rho_k_1 = st.session_state.vrho_k_1, alfa_1_1 = st.session_state.valfa_1_1, fi_1 = st.session_state.vfi_1,
                        rho_k_z = st.session_state.vrho_k_z, alfa_1_z = st.session_state.valfa_1_z, fi_z = st.session_state.vfi_z, 
                        D_sr_reg = st.session_state.vD_sr_reg, delta_D = st.session_state.vdelta_D, delta = st.session_state.vdelta,
                        n = st.session_state.vn)
                        st.session_state.bdown = bdown
                        Save_to_file_stage(breakdown_table(st.session_state.bdown[0], 'param_1'), name = 'Результаты расчета ЦВД', extension = '.xlsx')
                        Save_to_file_stage(breakdown_table(st.session_state.bdown[2], 'param_2'), name = 'Распределение параметров по ступеням', extension = '.xlsx')

        if select_1 ==  'Расчет ЦСД':
            with st.form(key = 'my_form_2'):
                    _aP_0_, _at_0_, _aP_2_z, = st.columns(3)
                    _aG_0, _aG_z, _aetta_oi = st.columns(3)
                    _arho_k_1, _aalfa_1_1, _afi_1, = st.columns(3)
                    _arho_k_z, _aalfa_1_z, _afi_z = st.columns(3)
                    _aD_k_2_1, _aD_k_2_z, _adelta, = st.columns(3)
                    _an, _ = st.columns(2)
                    st.session_state.method_form_s = 0.0
                    st.session_state.ai = 1

                    P_0_ = count(counter = 'aP_0_', column = _aP_0_, name = "P_0, МПа", format = "%f", ke = 'count18')
                    t_0_ = count(counter = 'at_0_', column = _at_0_, name = "t_0, ℃", format = "%f", ke = 'count19')
                    P_2_z = count(counter = 'aP_2_z', column = _aP_2_z, name = "P_2z, МПа", format = "%f", ke = 'count20')
                    G_0 = count(counter = 'aG_0', column = _aG_0, name = "G_0, кг/с", format = "%f", ke = 'count21')
                    G_z = count(counter = 'aG_z', column = _aG_z, name = "G_z, кг/с", format = "%f", ke = 'count22')
                    rho_k_1 = count(counter = 'arho_k_1', column = _arho_k_1, name = "rho_k_1, -", format = "%f", ke = 'count24')
                    alfa_1_1 = count(counter = 'aalfa_1_1', column = _aalfa_1_1, name = "alfa_1_1, град", format = "%f", ke = 'count25')
                    fi_1 = count(counter = 'afi_1', column = _afi_1, name = "fi_1, -", format = "%f", ke = 'count26')
                    rho_k_z = count(counter = 'arho_k_z', column = _arho_k_z, name = "rho_k_z, -", format = "%f", ke = 'count27')
                    alfa_1_z = count(counter = 'aalfa_1_z', column = _aalfa_1_z, name = "alfa_1_z, град", format = "%f", ke = 'count28')
                    fi_z = count(counter = 'afi_z', column = _afi_z, name = "fi_z, -", format = "%f", ke = 'count29')
                    D_k_2_1 = count(counter = 'aD_k_2_1', column = _aD_k_2_1, name = "D_k_2_1, м", format = "%f", ke = 'count30')
                    D_k_2_z = count(counter = 'aD_k_2_z', column = _aD_k_2_z, name = "D_k_2_z, м", format = "%f", ke = 'count31')
                    delta = count(counter = 'adelta', column = _adelta, name = "delta, м", format = "%f", ke = 'count32')
                    etta_oi = count(counter = 'aetta_oi', column = _aetta_oi, name = "etta_oi, -", format = "%f", ke = 'count33')
                    n = count(counter = 'an', column = _an, name = "n, Гц", format = "%g", ke = 'count34')
                    
                    # st.write('P_0_ = ',P_0_, '; ','t_0_ = ',t_0_, '; ','P_2_z = ',P_2_z)
                    # st.write('G_0 = ',G_0, '; ','G_z = ',G_z, '; ','mu_1 = ',mu_1)
                    # st.write('rho_k_1 = ',rho_k_1, '; ','alfa_1_1 = ',alfa_1_1, '; ','fi_1 = ',fi_1)
                    # st.write('rho_k_z = ',rho_k_z, '; ','alfa_1_z = ',alfa_1_z, '; ','fi_z = ',fi_z)
                    # st.write('D_k_2_1 = ',D_k_2_1, '; ','D_k_2_z = ',D_k_2_z, '; ','delta = ',delta)
                    # st.write('etta_oi = ',etta_oi, '; ','n = ',n)
                    
                    form, _ai, = st.columns(2)
                    with form:
                        form = st.selectbox('Форма проточной части', options = ('Форма 1', 'Форма 2', 'Форма 3', 'Форма 4'))
                        i_ = count(counter = 'ai', column = _ai, name = "i, -", format = "%g", ke = 'count35')

                    if form == 'Форма 1':
                        st.session_state.method_form_s = 'form_1'
                    if form == 'Форма 2':
                        st.session_state.method_form_s = 'form_2'
                    if form == 'Форма 3':
                        st.session_state.method_form_s = 'form_3'
                    if form == 'Форма 4':
                        st.session_state.method_form_s = 'form_4'

                    if st.form_submit_button('Расчет '):
                        bdown = breakdown_CSP(P_0_ = st.session_state.aP_0_, t_0_ = st.session_state.at_0_, P_2_z =  st.session_state.aP_2_z, 
                        G_0 = st.session_state.aG_0, G_z = st.session_state.aG_z, etta_oi = st.session_state.aetta_oi,
                        rho_k_1 = st.session_state.arho_k_1, alfa_1_1 = st.session_state.aalfa_1_1, fi_1 = st.session_state.afi_1,
                        rho_k_z = st.session_state.arho_k_z, alfa_1_z = st.session_state.aalfa_1_z, fi_z = st.session_state.afi_z, 
                        D_k_2_1 = st.session_state.aD_k_2_1, D_k_2_z = st.session_state.aD_k_2_z, delta = st.session_state.adelta,
                        n = st.session_state.an, method_1= st.session_state.method_form_s, i = st.session_state.ai)
                        st.session_state.bdown = bdown
                        st.header('Результаты расчета ЦСД')
                        st.table(breakdown_table(st.session_state.bdown[0], 'param_1'))
                        st.header('Распределение параметров по ступеням')
                        st.table(breakdown_table(st.session_state.bdown[2], 'param_2'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][4], k =  st.session_state.bdown[2][3], z = st.session_state.bdown[2][2],  method = 'Form_diametr'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][3], method = 'Diametr'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][1], method = 'Height'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][5], method = 'Tetta'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][7], method = 'Reactivity'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][8], method = 'X'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][11], method = 'Heat_difference'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][12], method = 'Heat_difference_new'))
                    
                    if st.form_submit_button('Сохранить в Excel'):
                        bdown = breakdown_CSP(P_0_ = st.session_state.aP_0_, t_0_ = st.session_state.at_0_, P_2_z =  st.session_state.aP_2_z, 
                        G_0 = st.session_state.aG_0, G_z = st.session_state.aG_z, etta_oi = st.session_state.aetta_oi,
                        rho_k_1 = st.session_state.arho_k_1, alfa_1_1 = st.session_state.aalfa_1_1, fi_1 = st.session_state.afi_1,
                        rho_k_z = st.session_state.arho_k_z, alfa_1_z = st.session_state.aalfa_1_z, fi_z = st.session_state.afi_z, 
                        D_k_2_1 = st.session_state.aD_k_2_1, D_k_2_z = st.session_state.aD_k_2_z, delta = st.session_state.adelta,
                        n = st.session_state.an, method_1= st.session_state.method_form_s, i = st.session_state.ai)
                        st.session_state.bdown = bdown
                        Save_to_file_stage(breakdown_table(st.session_state.bdown[0], 'param_1'), name = 'Результаты расчета ЦСД', extension = '.xlsx')
                        Save_to_file_stage(breakdown_table(st.session_state.bdown[2], 'param_2'), name = 'Распределение параметров по ступеням', extension = '.xlsx')

        if select_1 ==  'Расчет ЦНД':
            with st.form(key = 'my_form_3'):
                    _lP_0_, _lt_0_, _lP_2_z, = st.columns(3)
                    _lG_0, _lG_z, _letta_oi = st.columns(3)
                    _lrho_k_1, _lalfa_1_1, _lfi_1, = st.columns(3)
                    _lrho_k_z, _lalfa_1_z, _lfi_z = st.columns(3)
                    _lD_k_2_1, _lD_k_2_z, _ldelta, = st.columns(3)
                    _ln, _ = st.columns(2)
                    st.session_state.method_form_l = 0.0
                    st.session_state.li = 1

                    P_0_ = count(counter = 'lP_0_', column = _lP_0_, name = "P_0, МПа", format = "%f", ke = 'count36')
                    t_0_ = count(counter = 'lt_0_', column = _lt_0_, name = "t_0, ℃", format = "%f", ke = 'count37')
                    P_2_z = count(counter = 'lP_2_z', column = _lP_2_z, name = "P_2z, МПа", format = "%f", ke = 'count38')
                    G_0 = count(counter = 'lG_0', column = _lG_0, name = "G_0, кг/с", format = "%f", ke = 'count39')
                    G_z = count(counter = 'lG_z', column = _lG_z, name = "G_z, кг/с", format = "%f", ke = 'count40')
                    rho_k_1 = count(counter = 'lrho_k_1', column = _lrho_k_1, name = "rho_k_1, -", format = "%f", ke = 'count41')
                    alfa_1_1 = count(counter = 'lalfa_1_1', column = _lalfa_1_1, name = "alfa_1_1, град", format = "%f", ke = 'count42')
                    fi_1 = count(counter = 'lfi_1', column = _lfi_1, name = "fi_1, -", format = "%f", ke = 'count43')
                    rho_k_z = count(counter = 'lrho_k_z', column = _lrho_k_z, name = "rho_k_z, -", format = "%f", ke = 'count44')
                    alfa_1_z = count(counter = 'lalfa_1_z', column = _lalfa_1_z, name = "alfa_1_z, град", format = "%f", ke = 'count45')
                    fi_z = count(counter = 'lfi_z', column = _lfi_z, name = "fi_z, -", format = "%f", ke = 'count46')
                    D_k_2_1 = count(counter = 'lD_k_2_1', column = _lD_k_2_1, name = "D_k_2_1, м", format = "%f", ke = 'count47')
                    D_k_2_z = count(counter = 'lD_k_2_z', column = _lD_k_2_z, name = "D_k_2_z, м", format = "%f", ke = 'count48')
                    delta = count(counter = 'ldelta', column = _ldelta, name = "delta, м", format = "%f", ke = 'count49')
                    etta_oi = count(counter = 'letta_oi', column = _letta_oi, name = "etta_oi, -", format = "%f", ke = 'count50')
                    n = count(counter = 'ln', column = _ln, name = "n, Гц", format = "%g", ke = 'count51')
                    
                    # st.write('P_0_ = ',P_0_, '; ','t_0_ = ',t_0_, '; ','P_2_z = ',P_2_z)
                    # st.write('G_0 = ',G_0, '; ','G_z = ',G_z, ';)
                    # st.write('rho_k_1 = ',rho_k_1, '; ','alfa_1_1 = ',alfa_1_1, '; ','fi_1 = ',fi_1)
                    # st.write('rho_k_z = ',rho_k_z, '; ','alfa_1_z = ',alfa_1_z, '; ','fi_z = ',fi_z)
                    # st.write('D_k_2_1 = ',D_k_2_1, '; ','D_k_2_z = ',D_k_2_z, '; ','delta = ',delta)
                    # st.write('etta_oi = ',etta_oi, '; ','n = ',n)
                    
                    form, _li, = st.columns(2)
                    with form:
                        form = st.selectbox('Форма проточной части', options = ('Форма 1', 'Форма 2', 'Форма 3', 'Форма 4'))
                        i_ = count(counter = 'li', column = _li, name = "i, -", format = "%g", ke = 'count52')

                    if form == 'Форма 1':
                        st.session_state.method_form_l = 'form_1'
                    if form == 'Форма 2':
                        st.session_state.method_form_l = 'form_2'
                    if form == 'Форма 3':
                        st.session_state.method_form_l = 'form_3'
                    if form == 'Форма 4':
                        st.session_state.method_form_l = 'form_4'

                    if st.form_submit_button('Расчет '):
                        bdown = breakdown_CSP(P_0_ = st.session_state.lP_0_, t_0_ = st.session_state.lt_0_, P_2_z =  st.session_state.lP_2_z, 
                        G_0 = st.session_state.lG_0, G_z = st.session_state.lG_z, etta_oi = st.session_state.letta_oi,
                        rho_k_1 = st.session_state.lrho_k_1, alfa_1_1 = st.session_state.lalfa_1_1, fi_1 = st.session_state.lfi_1,
                        rho_k_z = st.session_state.lrho_k_z, alfa_1_z = st.session_state.lalfa_1_z, fi_z = st.session_state.lfi_z, 
                        D_k_2_1 = st.session_state.lD_k_2_1, D_k_2_z = st.session_state.lD_k_2_z, delta = st.session_state.ldelta,
                        n = st.session_state.ln, method_1= st.session_state.method_form_l, i = st.session_state.li)
                        st.session_state.bdown = bdown
                        
                        st.header('Результаты расчета ЦНД')
                        st.table(breakdown_table(st.session_state.bdown[0], 'param_1'))
                        st.header('Распределение параметров по ступеням')
                        st.table(breakdown_table(st.session_state.bdown[2], 'param_2'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][4], k =  st.session_state.bdown[2][3], z = st.session_state.bdown[2][2],  method = 'Form_diametr'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][3], method = 'Diametr'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][1], method = 'Height'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][5], method = 'Tetta'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][7], method = 'Reactivity'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][8], method = 'X'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][11], method = 'Heat_difference'))
                        st.pyplot(breakdown_plot(x = st.session_state.bdown[2][0], y = st.session_state.bdown[2][12], method = 'Heat_difference_new'))
                    
                    if st.form_submit_button('Сохранить в Excel'):
                        bdown = breakdown_CSP(P_0_ = st.session_state.lP_0_, t_0_ = st.session_state.lt_0_, P_2_z =  st.session_state.lP_2_z, 
                        G_0 = st.session_state.lG_0, G_z = st.session_state.lG_z, etta_oi = st.session_state.letta_oi,
                        rho_k_1 = st.session_state.lrho_k_1, alfa_1_1 = st.session_state.lalfa_1_1, fi_1 = st.session_state.lfi_1,
                        rho_k_z = st.session_state.lrho_k_z, alfa_1_z = st.session_state.lalfa_1_z, fi_z = st.session_state.lfi_z, 
                        D_k_2_1 = st.session_state.lD_k_2_1, D_k_2_z = st.session_state.lD_k_2_z, delta = st.session_state.ldelta,
                        n = st.session_state.ln, method_1= st.session_state.method_form_l, i = st.session_state.li)
                        st.session_state.bdown = bdown
                        Save_to_file_stage(breakdown_table(st.session_state.bdown[0], 'param_1'), name = 'Результаты расчета ЦСД', extension = '.xlsx')
                        Save_to_file_stage(breakdown_table(st.session_state.bdown[2], 'param_2'), name = 'Распределение параметров по ступеням', extension = '.xlsx')

    if panel_1 == "2. - Расчет геометрии ПТУ":
        st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет геометрии ПТУ</ins> </h1>", unsafe_allow_html=True)
        
        if st.session_state.bdown == None:
            st.header('Необходимо выполнить: "1. - Расчет цилиндров ПТУ"')
        else:
            st.header('Ввод исходных данных')
            with st.form(key = 'my_form_4'):
                _K_s, _K_r, _a_c, = st.columns(3)

                K_s = count(counter = 'K_s', column = _K_s, name = "Коэффициент ширины сопловой, -", format = "%f", ke = 'count53')
                K_r = count(counter = 'K_r', column = _K_r, name = "Коэффициент ширины рабочей, -", format = "%f", ke = 'count54')
                a_c = count(counter = 'r_c', column = _a_c, name = "Коэффициент осевого зазор, -", format = "%f", ke = 'count55')
                
                if st.form_submit_button('Расчет  '): 
                    geom = geometry(br = st.session_state.bdown, K_s = st.session_state.K_s, K_r = st.session_state.K_r, axial_clearance = st.session_state.r_c)
                    st.session_state.geom = geom
                    st.pyplot(flowpath_plot(st.session_state.geom))
                    st.header('Результаты расчета геометрии проточной части')
                    st.table(geometry_table(st.session_state.geom[2]))

                if st.form_submit_button('Сохранить в Excel'):
                    geom = geometry(br = st.session_state.bdown, K_s = st.session_state.K_s, K_r = st.session_state.K_r, axial_clearance = st.session_state.r_c)
                    st.session_state.geom = geom
                    Save_to_file_stage(geometry_table(st.session_state.geom[2]), name = 'Результаты расчета геометрии проточной части', extension = '.xlsx')

    if panel_1 == "3. - Расчет параметров по ступеням":
        st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет параметров по ступеням</ins> </h1>", unsafe_allow_html=True)
        if st.session_state.bdown == None:
            st.header('Необходимо выполнить: "1. - Расчет цилиндров ПТУ"')
        else:
            with st.form(key = 'my_form_5'):
                if st.form_submit_button('Расчет'):
                    heat = heattransfer(br = st.session_state.bdown)
                    st.session_state.heat = heat
                    st.session_state.number_of_steps_ = st.session_state.heat[0]['number_of_steps']
                    st.pyplot(hs_plot(point_0_t_i_ = st.session_state.heat[2][0], point_2_t_i_ = st.session_state.heat[2][1], point_0_i = st.session_state.heat[2][2], point_2_i = st.session_state.heat[2][3], point_0_i_ = st.session_state.heat[2][4], point_2_i_ = st.session_state.heat[2][5]))
                    st.header('Результаты расчета параметров по ступеням')
                    st.table(heattransfer_table(st.session_state.heat[1]))          
                
                if st.form_submit_button('Сохранить в Excel'):
                    heat = heattransfer(br = st.session_state.bdown)
                    st.session_state.heat = heat
                    Save_to_file_stage(heattransfer_table(st.session_state.heat[1]), name = 'Результаты расчета параметров по ступеням', extension = '.xlsx')

if panel_global == "II. - Этап расчета ступеней ПТУ":
    if st.session_state.bdown == None:
        st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет ступени ПТУ на среднем диаметре</ins></h1>", unsafe_allow_html=True)
        st.header('Необходимо выполнить: "I. - Этап расчета ПТУ"; "3 - Расчет параметров по ступеням"')
    else:
        stage_str = [f'Ступень №{i+1}' for i in range(int(st.session_state.number_of_steps_))]
        stage_str.append('Параметры расчета по ступеням')
        panel_2 = st.sidebar.radio('Этапы расчета ступени:', stage_str)
        
        num = [f'Ступень №{i+1}' for i in range(int(st.session_state.number_of_steps_))]
        for i in range(int(st.session_state.number_of_steps_)):
            if panel_2 == num[i]:
                st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет ступени ПТУ на среднем диаметре</ins></h1>", unsafe_allow_html=True)
                st.session_state.stages = None
                st.header(f'Ввод исходных данных ступени №{i+1}')
                with st.form(key = 'my_form_6'): 
                    if i == 0:
                        st.session_state.P0 = []
                        st.session_state.t0 = []
                        st.session_state.x0 = []
                        st.session_state.C0 = []
                        st.session_state.G0 = []
                        st.session_state.alpha2 = []
                        st.session_state.paramstage = []
                    st.session_state.value1_sopl_, st.session_state.value1_rab_ = 0.0, 0.0
                    st.session_state.value2_sopl_, st.session_state.value2_rab_ = 0.0, 0.0
                    st.session_state.value3_sopl_, st.session_state.value3_rab_ = 0.0, 0.0

                    st.header('Параметры сопловой:') 
                    param_2 = countSRVVV(i, 'sopl')   
                    st.session_state.value1_sopl_ = param_2[0]
                    st.session_state.value2_sopl_ = param_2[1]
                    st.session_state.value3_sopl_ = param_2[2]
                    losses_sopl = st.selectbox('Модель потерь для сопловой решетки', options = ('Soderberg losses', 'Ainley and Mathieson losses', 'Craig and Cox losses', 'CIAM losses', 'Denton losses'))
                            
                    if st.form_submit_button('Выбрать '):
                        st.session_state.method_losses_sopl = True
                        st.session_state.coef_sopl_ = 0.0
                        st.session_state.B_sopl_ = 0.0
                        st.session_state.ks_sopl_ = 0.0
                        st.session_state.sorU_sopl_ = 0.0

                    if losses_sopl == 'Soderberg losses':
                        st.session_state.method_losses_sopl = 'SDB'

                    if losses_sopl == 'Ainley and Mathieson losses':
                        param_3 = countLosses(i,'ANM', 'sopl')   
                        st.session_state.coef_sopl_ = param_3[0]
                        st.session_state.B_sopl_ = param_3[1]
                        st.session_state.method_losses_sopl = 'ANM'

                    if losses_sopl == 'Craig and Cox losses':
                        param_3 = countLosses(i,'CAC', 'sopl')   
                        st.session_state.ks_sopl_ = param_3[0]
                        st.session_state.method_losses_sopl = 'CAC'

                    if losses_sopl == 'CIAM losses':
                        st.session_state.method_losses_sopl = 'CIAM'

                    if losses_sopl == 'Denton losses':
                        param_3 = countLosses(i,'DN', 'sopl')   
                        st.session_state.sorU_sopl_ = param_3[0]
                        st.session_state.method_losses_sopl = 'DN'

                    st.header('Параметры рабочей:') 
                    param_4 = countSRVVV(i, 'rab') 
                    st.session_state.value1_rab_ = param_4[0]
                    st.session_state.value2_rab_ = param_4[1]
                    st.session_state.value3_rab_ = param_4[2]
                    losses_rab = st.selectbox('Модель потерь для рабочей решетки', options = ('Soderberg losses', 'Ainley and Mathieson losses', 'Craig and Cox losses', 'CIAM losses', 'Denton losses'))
                    
                    if st.form_submit_button('Выбрать  '):
                        st.session_state.method_losses_rab = True
                        st.session_state.coef_rab_ = 0.0
                        st.session_state.B_rab_ = 0.0
                        st.session_state.ks_rab_ = 0.0
                        st.session_state.sorU_rab_ = 0.0  

                    if losses_rab == 'Soderberg losses':
                        st.session_state.method_losses_rab = 'SDB' 

                    if losses_rab == 'Ainley and Mathieson losses':
                        param_4 = countLosses(i,'ANM', 'rab')   
                        st.session_state.coef_rab_ = param_4[0]
                        st.session_state.B_rab_ = param_4[1]
                        st.session_state.method_losses_rab = 'ANM'

                    if losses_rab == 'Craig and Cox losses':
                        param_4 = countLosses(i,'CAC', 'rab')   
                        st.session_state.ks_rab_ = param_4[0]
                        st.session_state.method_losses_rab = 'CAC'
                            
                    if losses_rab == 'CIAM losses':
                        st.session_state.method_losses_rab = 'CIAM'

                    if losses_rab == 'Denton losses':
                        param_4 = countLosses(i,'DN', 'rab')   
                        st.session_state.sorU_rab_ = param_4[0]
                        st.session_state.method_losses_rab = 'DN'
                    
                    if st.form_submit_button('Расчет'):
                        if i == 0:
                            stages = stage(br = st.session_state.bdown, geom = st.session_state.geom, heat = st.session_state.heat, 
                            P_0 = st.session_state.heat[0]['P_0_i'][i], t_0 = st.session_state.heat[i]['t_0_i'][i], x_0 = st.session_state.heat[i]['x_0_i'][i],
                            C_0 = st.session_state.heat[0]['C_0_i'][i], G_0 = st.session_state.bdown[0]['G_0'], alpha_0 = 90, 
                            value1_sopl = st.session_state.value1_sopl_, value1_rab = st.session_state.value1_rab_ , value2_sopl = st.session_state.value2_sopl_, value2_rab = st.session_state.value2_rab_ , value3_sopl = st.session_state.value3_sopl_, value3_rab = st.session_state.value3_rab_ ,
                            coef_sopl = st.session_state.coef_sopl_, coef_rab = st.session_state.coef_rab_, B_sopl = st.session_state.B_sopl_, B_rab = st.session_state.B_rab_, SorU_sopl = st.session_state.sorU_sopl_, SorU_rab = st.session_state.sorU_rab_, 
                            ks_sopl = st.session_state.ks_sopl_, ks_rab = st.session_state.ks_rab_, method_losses_sopl = st.session_state.method_losses_sopl, method_losses_rab = st.session_state.method_losses_rab, num = i)
                            st.session_state.stages = stages 
                            
                            st.session_state.P0 = [] 
                            st.session_state.P0.insert(0, st.session_state.stages[5])  
  
                            st.session_state.t0 = [] 
                            st.session_state.t0.insert(0, st.session_state.stages[6])  

                            st.session_state.x0 = [] 
                            st.session_state.x0.insert(0, st.session_state.stages[7])

                            st.session_state.C0 = [] 
                            st.session_state.C0.insert(0, st.session_state.stages[8])  

                            st.session_state.G0 = [] 
                            st.session_state.G0.insert(0, st.session_state.stages[9]) 

                            st.session_state.alpha2 = [] 
                            st.session_state.alpha2.insert(0, st.session_state.stages[10])

                            st.session_state.paramstage = []
                            st.session_state.paramstage.insert(0, st.session_state.stages[0])
                            # st.write(st.session_state.P0)
                            # st.write(st.session_state.t0)
                            # st.write(st.session_state.C0) 
                            # st.write(st.session_state.G0)  
                            # st.write(st.session_state.alpha2)
                            # st.write(st.session_state.x0)
                            # st.write(st.session_state.paramstage)  
                        else:
                            stages = stage(br = st.session_state.bdown, geom = st.session_state.geom, heat = st.session_state.heat, 
                            P_0 = st.session_state.P0[i-1], t_0 = st.session_state.t0[i-1], x_0 = st.session_state.x0[i-1],
                            C_0 = st.session_state.C0[i-1], G_0 = st.session_state.G0[i-1], alpha_0 = st.session_state.alpha2[i-1], 
                            value1_sopl = st.session_state.value1_sopl_, value1_rab = st.session_state.value1_rab_ , value2_sopl = st.session_state.value2_sopl_, value2_rab = st.session_state.value2_rab_ , value3_sopl = st.session_state.value3_sopl_, value3_rab = st.session_state.value3_rab_ ,
                            coef_sopl = st.session_state.coef_sopl_, coef_rab = st.session_state.coef_rab_, B_sopl = st.session_state.B_sopl_, B_rab = st.session_state.B_rab_, SorU_sopl = st.session_state.sorU_sopl_, SorU_rab = st.session_state.sorU_rab_, 
                            ks_sopl = st.session_state.ks_sopl_, ks_rab = st.session_state.ks_rab_, method_losses_sopl = st.session_state.method_losses_sopl, method_losses_rab = st.session_state.method_losses_rab, num = i)
                            st.session_state.stages = stages 

                            if len(st.session_state.P0) == i: 
                                st.session_state.P0.insert(i, st.session_state.stages[5])
                                st.session_state.t0.insert(i, st.session_state.stages[6])
                                st.session_state.x0.insert(i, st.session_state.stages[7])
                                st.session_state.C0.insert(i, st.session_state.stages[8])
                                st.session_state.G0.insert(i, st.session_state.stages[9])
                                st.session_state.alpha2.insert(i, st.session_state.stages[10])
                                st.session_state.paramstage.insert(i, st.session_state.stages[0])
                            else: 
                                st.session_state.P0[i] = st.session_state.stages[5]
                                st.session_state.t0[i] = st.session_state.stages[6]
                                st.session_state.x0[i] = st.session_state.stages[7]
                                st.session_state.C0[i] = st.session_state.stages[8]
                                st.session_state.G0[i] = st.session_state.stages[9]
                                st.session_state.alpha2[i] = st.session_state.stages[10]
                                st.session_state.paramstage[i] = st.session_state.stages[0]          
                            # st.write(st.session_state.P0) 
                            # st.write(st.session_state.t0) 
                            # st.write(st.session_state.C0) 
                            # st.write(st.session_state.G0)
                            # st.write(st.session_state.alpha2)
                            # st.write(st.session_state.x0)       
                            # st.write(st.session_state.paramstage)
                        st.pyplot(hs_stage_plot(point_0_ = st.session_state.stages[1][0], point_0 = st.session_state.stages[1][1], point_1t = st.session_state.stages[1][2], point_1 = st.session_state.stages[1][3], 
                        point_1w = st.session_state.stages[1][4], point_2t_ = st.session_state.stages[1][5], point_2t = st.session_state.stages[1][6], point_2 = st.session_state.stages[1][7], point_2w = st.session_state.stages[1][8], 
                        delta_Htr = st.session_state.stages[1][9], delta_Hlake = st.session_state.stages[1][10], delta_Hvet = st.session_state.stages[1][11], Delta_Hvs = st.session_state.stages[1][12], kappa_vs = st.session_state.stages[1][13], num = i)) 
                        
                        st.pyplot(velocity_triangle_plot(C_1 = st.session_state.stages[2][0], W_1 = st.session_state.stages[2][1], U_1 = st.session_state.stages[2][2], alpha_1 = st.session_state.stages[2][3], 
                        betta_1 = st.session_state.stages[2][4], C_2 = st.session_state.stages[2][5], W_2 = st.session_state.stages[2][6], U_2 = st.session_state.stages[2][7], alpha_2 = st.session_state.stages[2][8], betta_2 = st.session_state.stages[2][9], num = i))
                
                        st.header('Термодинамические и кинематические параметры')
                        st.table(stage_table(st.session_state.stages[0], method = 'parameters'))
                        st.header('Геометрические параметры профиля сопловой')
                        st.table(stage_table(st.session_state.stages[3], method = "profile sopl"))
                        st.header('Геометрические параметры профиля рабочей')
                        st.table(stage_table(st.session_state.stages[4], method = "profile rab"))

                        Save_to_file_stage(stage_table(st.session_state.stages[0], method = 'parameters'), name = f'Ступень №{i+1} (Параметры на среднем диаметре)', extension = '.xlsx')
                        Save_to_file_stage(stage_table(st.session_state.stages[3], method = "profile sopl"), name = f'Соплова геом. №{i+1}', extension = '.xlsx')
                        Save_to_file_stage(stage_table(st.session_state.stages[4], method = "profile rab"), name = f'Рабочая геом. №{i+1}', extension = '.xlsx')

        if panel_2 == 'Параметры расчета по ступеням':
            st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Результаты расчета на среднем диаметре</ins></h1>", unsafe_allow_html=True)
            with st.form(key = 'my_form_7'): 
                st.session_state.stage_dict = {}
                for n in st.session_state.paramstage[0].keys():
                    st.session_state.stage_dict[n] = list(d[n] for d in st.session_state.paramstage)               
                st.session_state.stage_list = list(st.session_state.stage_dict.values()) 
                st.header('Результаты расчета параметров по ступеням') 
                st.table(stageTable(st.session_state.stage_list))
                if st.form_submit_button('Сохранить в Excel'):
                    Save_to_file_stage(stageTable(st.session_state.stage_list), name = 'Итоговые результаты расчета параметров по ступеням', extension = '.xlsx')

if panel_global == "III. - Этап расчета по сечениям":
        st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет ступени ПТУ по сечениям</ins></h1>", unsafe_allow_html=True)
        if st.session_state.bdown == None:
            st.header('Необходимо выполнить: "II. - Этап расчета ступеней ПТУ"; "3 - Параметры расчета по ступеням"')
        else:
            panel_3 = st.sidebar.radio('Этапы расчета ступени:', [f'Ступень №{i+1}' for i in range(int(st.session_state.number_of_steps_))])
            num_1 = [f'Ступень №{i+1}' for i in range(int(st.session_state.number_of_steps_))]
        
            for i in range(int(st.session_state.number_of_steps_)):
                    if panel_3 == num_1[i]:

                        st.session_state.section = None
                        st.session_state.method_section = True
                        st.session_state.value_num_  = 0.0
                        st.header(f'Ввод исходных данных ступени №{i+1}')

                        with st.form(key = 'my_form_8'): 
                            select_2 = st.radio('Выбор закона закрутки:', ['Закон постоянства угла выхода: 𝛼1(𝑟) = 𝑐𝑜𝑛𝑠𝑡', 'Обратный закон закрутки: 𝑟 ∙ 𝑡𝑔(𝛼1) = 𝑐𝑜𝑛𝑠𝑡', 'Закон постоянства циркуляции: 𝐶1𝑢 ∙ 𝑟𝜑2 = 𝑐𝑜𝑛𝑠𝑡'])  
                            if i == 0:
                                st.session_state.data = []
                            if select_2 == 'Обратный закон закрутки: 𝑟 ∙ 𝑡𝑔(𝛼1) = 𝑐𝑜𝑛𝑠𝑡':
                                param_5 = countNum(i)
                                st.session_state.value_num_ = param_5
                                st.session_state.method_section = 'rtgconst' 

                            if select_2 == 'Закон постоянства циркуляции: 𝐶1𝑢 ∙ 𝑟𝜑2 = 𝑐𝑜𝑛𝑠𝑡':
                                param_5 = countNum(i)
                                st.session_state.value_num_ = param_5
                                st.session_state.method_section = 'C1uconst'                   

                            if select_2 == 'Закон постоянства угла выхода: 𝛼1(𝑟) = 𝑐𝑜𝑛𝑠𝑡':
                                param_5 = countNum(i)
                                st.session_state.value_num_ = param_5
                                st.session_state.method_section = 'alpha1const'     

                            if st.form_submit_button('Расчет'):                        
                                section = spin_laws_stage(stg = st.session_state.stage_dict, num = i, sect = int(st.session_state.value_num_), method = st.session_state.method_section)
                                st.session_state.section = section
                                st.header(f'Результаты расчета параметров по сечениям ступени №{i+1}')
                                st.table(sectionTable(st.session_state.section[1]))
                                
                                st.pyplot(velocity_triangle_i(C_1_i = st.session_state.section[2][0], W_1_i = st.session_state.section[2][1], U_1_i = st.session_state.section[2][2], alpha_1_i = st.session_state.section[2][3], betta_1_i = st.session_state.section[2][4],
                                C_2_i = st.session_state.section[2][5], W_2_i = st.session_state.section[2][6], U_2_i = st.session_state.section[2][7], alpha_2_i = st.session_state.section[2][8], betta_2_i = st.session_state.section[2][9], num = i))
                                
                                st.pyplot(parametrs(alpha_1_i = st.session_state.section[2][3], betta_1_i = st.session_state.section[2][4], 
                                                    alpha_2_i = st.session_state.section[2][8], betta_2_i = st.session_state.section[2][9], 
                                                    C_1_i = st.session_state.section[2][0], C_2_i = st.session_state.section[2][5], 
                                                    W_1_i = st.session_state.section[2][1], W_2_i = st.session_state.section[2][6], 
                                                    M_1c_i = st.session_state.section[2][10], M_2c_i = st.session_state.section[2][11], 
                                                    M_1w_i = st.session_state.section[2][12], M_2w_i = st.session_state.section[2][13], 
                                                    fi_i = st.session_state.section[2][14], psi_i = st.session_state.section[2][15], 
                                                    num = i, sect = st.session_state.section[2][16], method = 'angle'))

                                st.pyplot(parametrs(alpha_1_i = st.session_state.section[2][3], betta_1_i = st.session_state.section[2][4], 
                                                    alpha_2_i = st.session_state.section[2][8], betta_2_i = st.session_state.section[2][9], 
                                                    C_1_i = st.session_state.section[2][0], C_2_i = st.session_state.section[2][5], 
                                                    W_1_i = st.session_state.section[2][1], W_2_i = st.session_state.section[2][6], 
                                                    M_1c_i = st.session_state.section[2][10], M_2c_i = st.session_state.section[2][11], 
                                                    M_1w_i = st.session_state.section[2][12], M_2w_i = st.session_state.section[2][13], 
                                                    fi_i = st.session_state.section[2][14], psi_i = st.session_state.section[2][15], 
                                                    num = i, sect = st.session_state.section[2][16], method = 'velocity'))

                                st.pyplot(parametrs(alpha_1_i = st.session_state.section[2][3], betta_1_i = st.session_state.section[2][4], 
                                                    alpha_2_i = st.session_state.section[2][8], betta_2_i = st.session_state.section[2][9], 
                                                    C_1_i = st.session_state.section[2][0], C_2_i = st.session_state.section[2][5], 
                                                    W_1_i = st.session_state.section[2][1], W_2_i = st.session_state.section[2][6], 
                                                    M_1c_i = st.session_state.section[2][10], M_2c_i = st.session_state.section[2][11], 
                                                    M_1w_i = st.session_state.section[2][12], M_2w_i = st.session_state.section[2][13], 
                                                    fi_i = st.session_state.section[2][14], psi_i = st.session_state.section[2][15], 
                                                    num = i, sect = st.session_state.section[2][16], method = 'Mach'))

                                st.pyplot(parametrs(alpha_1_i = st.session_state.section[2][3], betta_1_i = st.session_state.section[2][4], 
                                                    alpha_2_i = st.session_state.section[2][8], betta_2_i = st.session_state.section[2][9], 
                                                    C_1_i = st.session_state.section[2][0], C_2_i = st.session_state.section[2][5], 
                                                    W_1_i = st.session_state.section[2][1], W_2_i = st.session_state.section[2][6], 
                                                    M_1c_i = st.session_state.section[2][10], M_2c_i = st.session_state.section[2][11], 
                                                    M_1w_i = st.session_state.section[2][12], M_2w_i = st.session_state.section[2][13], 
                                                    fi_i = st.session_state.section[2][14], psi_i = st.session_state.section[2][15], 
                                                    num = i, sect = st.session_state.section[2][16], method = 'losses'))
                            
                            if st.form_submit_button('Сохранить в Excel'):
                                section = spin_laws_stage(stg = st.session_state.stage_dict, num = i, sect = int(st.session_state.value_num_), method = st.session_state.method_section)
                                st.session_state.section = section
                                Save_to_file_stage(sectionTable(st.session_state.section[1]), name = f'Ступень №{i+1} (Параметры по сечениям)', extension = '.xlsx')



