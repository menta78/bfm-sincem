wavelenghts_list=["0250", "0325", "0350", "0375" ,"0400", "0425", "0450", "0475", "0500",
                  "0525", "0550", "0575", "0600", "0625", "0650", "0675", "0700", "0725",
                  "0775", "0850", "0950", "1050", "1150", "1250", "1350", "1450", "1550",
                  "1650", "1750", "1900", "2200", "2900", "3700" ]
      
### code block1 to be included in gotm.yaml
for i,wl in enumerate(wavelenghts_list):
    c=i+1
    string1="      surf_direct_downward_irradiance_" + wl + "_nm:"
    string2="         method: file                  # method [constant, file=from file; default=constant]"
    string3="         constant_value: 0.0         # value to use throughout the simulation [m/s; default=0.0]"
    string4="         file: GOTM_OASIM/Ed_" + wl + ".txt               # path to file with time series [default=]"
    string5="         column: 1                      # index of column to read from [default=1]"
    print(string1)
    print(string2)
    print(string3)
    print(string4)
    print(string5)

for i,wl in enumerate(wavelenghts_list):
    c=i+1
    string1="      surf_diffuse_downward_irradiance_" + wl + "_nm:"
    string2="         method: file                  # method [constant, file=from file; default=constant]"
    string3="         constant_value: 0.0         # value to use throughout the simulation [m/s; default=0.0]"
    string4="         file: GOTM_OASIM/Es_" + wl + ".txt               # path to file with time series [default=]"
    string5="         column: 1                     # index of column to read from [default=1]"
    print(string1)
    print(string2)
    print(string3)
    print(string4)
    print(string5)

