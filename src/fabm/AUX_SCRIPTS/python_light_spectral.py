def code_type_d(s,e):
    string="      type (type_horizontal_dependency_id) ::  id_Ed_0_" + wavelenghts_list[s] 
    for wl in wavelenghts_list[s+1:e]:
        string = string + ", id_Ed_0_" + wl
    print(string)


def code_type_s(s,e):
    string="      type (type_horizontal_dependency_id) ::  id_Es_0_" + wavelenghts_list[s] 
    for wl in wavelenghts_list[s+1:e]:
        string = string + ", id_Es_0_" + wl
    print(string)

wavelenghts_list=["0250", "0325", "0350", "0375" ,"0400", "0425", "0450", "0475", "0500",
                  "0525", "0550", "0575", "0600", "0625", "0650", "0675", "0700", "0725",
                  "0775", "0850", "0950", "1050", "1150", "1250", "1350", "1450", "1550",
                  "1650", "1750", "1900", "2200", "2900", "3700" ]
      
s="      "
### code block1 to be included in light_spectral.F90

code_type_d(0,5)
code_type_d(5,10)
code_type_d(10,15)
code_type_d(15,20)
code_type_d(20,25)
code_type_d(25,30)
code_type_d(30,33)

code_type_s(0,5)
code_type_s(5,10)
code_type_s(10,15)
code_type_s(15,20)
code_type_s(20,25)
code_type_s(25,30)
code_type_s(30,33)


### code block2 to be included in light_spectral.F90
for wl in wavelenghts_list:
      string1=s + "call self%register_dependency(self%id_Ed_0_" + wl + ",type_surface_standard_variable(name='surf_direct_downward_irradiance_" + wl + "_nm'))"
      print(string1)
for wl in wavelenghts_list:
      string2=s + "call self%register_dependency(self%id_Es_0_" + wl + ",type_surface_standard_variable(name='surf_diffuse_downward_irradiance_" + wl + "_nm'))"
      print(string2)

### code block3

for i,wl in enumerate(wavelenghts_list):
      string=s+ "_GET_SURFACE_(self%id_Ed_0_"  + wl +",Ed_0(" + str(i+1) + "))"
      print(string)

for i,wl in enumerate(wavelenghts_list):
      string=s+ "_GET_SURFACE_(self%id_Es_0_"  + wl +",Es_0(" + str(i+1) + "))"
      print(string)

