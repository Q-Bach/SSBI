# my personal settings
set stick_radius, 0.12
bg white
set auto_zoom, off
set cartoon_oval_length , 0.8
set cartoon_oval_width , 0.2
set ray_shadows,off
set ray_opaque_background, on
# zeigt doppelbindungen etc. an.
set valence, 0
set label_size, 30

fetch 1IGT

hide all
util.cbc
show cartoon, 1IGT
center 1IGT
set_view (\
    -0.269014597,   -0.650409758,    0.710342228,\
     0.861209214,    0.167761117,    0.479763448,\
    -0.431207538,    0.740815938,    0.515011191,\
     0.000210449,    0.000116058, -353.556274414,\
    12.967004776,  -21.872188568,    2.305969954,\
   306.241516113,  400.821777344,  -20.000000000 )
ray 2000,1500
png 1IGT.png, dpi=300

# LC1
hide all
create lc1, /1IGT/A/A
show cartoon, lc1
set_view (\
    -0.108022839,    0.461575925,   -0.880488217,\
     0.822625995,   -0.455802947,   -0.339872092,\
    -0.558206975,   -0.761030972,   -0.330465555,\
     0.001492269,    0.000361994, -154.142379761,\
    -7.286284447,  -53.447113037,   46.311168671,\
   129.062301636,  179.275451660,  -20.000000000 )
ray 2000,1500
png 1IGT_lc1.png, dpi=300

# HC1
hide all
create hc1, /1IGT/B/B
show cartoon, hc1
set_view (\
    -0.067964368,    0.760909140,   -0.645276010,\
     0.521201849,   -0.524414420,   -0.673292220,\
    -0.850707352,   -0.382083178,   -0.360945553,\
     0.000778869,   -0.000546496, -277.346313477,\
    10.040174484,  -44.209877014,   11.612754822,\
   244.589645386,  310.093353271,  -20.000000000 )
ray 2000,1500
png 1IGT_hc1.png, dpi=300

######## H-bridges in beta sheet of lc1
hide all
show sticks, lc1
color atomic, lc1
set_view (\
    -0.317352057,    0.040797636,    0.947410524,\
    -0.163517430,   -0.986425996,   -0.012292501,\
     0.934073687,   -0.158819631,    0.319721133,\
     0.002208221,   -0.002559205,  -24.974336624,\
    -4.623897552,  -42.711582184,   34.149810791,\
    17.270296097,   32.623901367,  -20.000000000 )

dist h01, /lc1/A/A/LEU`179/H, /lc1/A/A/VAL`132/O
dist h02, /lc1/A/A/CYS`134/H, /lc1/A/A/SER`177/O
dist h03, /lc1/A/A/SER`177/H, /lc1/A/A/CYS`134/O
group 1IGT_hbonds, h01 h02 h03
color black, 1IGT_hbonds
label n. CA, resn+resi
ray 2000,1500
png 1IGT_hbonds.png, dpi=300

######## disulfide bonds
hide all
show sticks, 1IGT
color atomic, 1IGT
set_view (\
    -0.280063063,    0.142434180,    0.949333668,\
    -0.606819332,   -0.792508662,   -0.060106490,\
     0.743819296,   -0.592920721,    0.308394283,\
    -0.001005173,   -0.000406715,  -37.053325653,\
    12.686450958,   -8.802280426,   12.922340393,\
    26.503728867,   48.065174103,  -20.000000000 )

dist s01, /1IGT/D/D/CYS`242/SG, /1IGT/B/B/CYS`242/SG
dist s02, /1IGT/D/D/CYS`240/SG, /1IGT/B/B/CYS`240/SG
dist s03, /1IGT/D/D/CYS`237/SG, /1IGT/B/B/CYS`237/SG
group 1IGT_sbonds, s01 s02 s03
color black, 1IGT_sbonds
ray 2000,1500
png 1IGT_sbonds.png, dpi=300

####### show all disulfides
hide all
show sticks, 1IGT
color atomic, 1IGT
set_view (\
    -0.269014597,   -0.650409758,    0.710342228,\
     0.861209214,    0.167761117,    0.479763448,\
    -0.431207538,    0.740815938,    0.515011191,\
     0.000210449,    0.000116058, -353.556274414,\
    12.967004776,  -21.872188568,    2.305969954,\
   306.241516113,  400.821777344,  -20.000000000 )
create 1IGT_disulfides, elem s and bound_to elem s and /1IGT
show sphere, 1IGT_disulfides
list=[]
iterate (1IGT_disulfides), list.append((model, chain, resn, resi))
python
for i in list:
    print(i[0], i[1], i[2], i[3])
python end

#############################################################
fetch 5IRE
util.cbc

# Complex
hide all
show cartoon, 5IRE
set_view (\
    -0.607543826,   -0.478832811,   -0.633713424,\
     0.774262786,   -0.535013497,   -0.338041902,\
    -0.177178487,   -0.696035326,    0.695789754,\
    -0.001428753,   -0.002847299, -403.669311523,\
  -113.445358276, -106.566207886, -122.098312378,\
   356.412994385,  450.993194580,  -20.000000000 )
ray 2000,1500
png 5IRE.png, dpi=300

# e-protein
create e_protein, /5IRE/A/A
hide all
show cartoon, e_protein
set_view (\
    -0.260012358,   -0.521334648,   -0.812764645,\
     0.869438469,   -0.492578059,    0.037808333,\
    -0.420061260,   -0.696817398,    0.581350148,\
    -0.001605265,   -0.004315868, -262.444854736,\
  -138.412597656,  -86.101104736, -116.154922485,\
   232.053253174,  292.811431885,  -20.000000000 )
ray 2000,1500
png 5IRE_e_protein.png, dpi=300

# m-protein
create m_protein, /5IRE/B/B
hide all
show cartoon, m_protein
set_view (\
     0.010377208,   -0.758987844,    0.651009202,\
    -0.440051794,   -0.588081360,   -0.678602517,\
     0.897898734,   -0.279434532,   -0.340098798,\
     0.002381410,   -0.002887856, -113.228614807,\
  -131.151535034,  -95.802688599, -101.110374451,\
    92.800773621,  133.473526001,  -20.000000000 )
ray 2000,1500
png 5IRE_m_protein.png, dpi=300

##### h-bonds in m-protein
hide all
show sticks, m_protein
color atomic, m_protein
h_add m_protein
set_view (\
     0.562865317,   -0.097675741,    0.820739686,\
     0.700932145,   -0.469801486,   -0.536617219,\
     0.438003004,    0.877328396,   -0.195974082,\
     0.005553079,    0.001618097,  -23.843645096,\
  -128.496368408,  -92.760383606,  -94.092071533,\
    14.802749634,   32.934810638,  -20.000000000 )
dist h04, /m_protein/B/B/ILE`67/H05, /m_protein/B/B/TYR`63/O
dist h05, /m_protein/B/B/MET`66/H10, /m_protein/B/B/ILE`62/O
dist h06, /m_protein/B/B/ILE`70/H12, /m_protein/B/B/MET`66/O
group 5IRE_hbonds, h04 h05 h06
color black, 5IRE_hbonds
label n. CA, resn+resi
ray 2000,1500
png 5IRE_hbonds.png, dpi=300

##### disulfide bonds in e-protein
hide all
show sticks, e_protein
color atomic, e_protein
set_view (\
     0.406716168,   -0.901010573,   -0.150693163,\
    -0.674181938,   -0.184745207,   -0.715022326,\
     0.616423965,    0.392416000,   -0.682600856,\
     0.001555835,   -0.003770374,  -82.799598694,\
  -126.119132996, -123.464408875, -101.776153564,\
  -146.014022827,  280.369140625,  -20.000000000 )
dist s04, /e_protein/A/A/CYS`105/SG, /e_protein/A/A/CYS`74/SG
dist s05, /e_protein/A/A/CYS`116/SG, /e_protein/A/A/CYS`92/SG
#dist s06, /5IRE/C/C/CYS`105/SG, /5IRE/C/C/CYS`74/SG
group 5IRE_sbonds, s04, s05, s06
color black, 5IRE_sbonds

####### show all disulfides
hide all
show sticks, 5IRE
color atomic, 5IRE
set_view (\
    -0.607543826,   -0.478832811,   -0.633713424,\
     0.774262786,   -0.535013497,   -0.338041902,\
    -0.177178487,   -0.696035326,    0.695789754,\
    -0.001428753,   -0.002847299, -403.669311523,\
  -113.445358276, -106.566207886, -122.098312378,\
   356.412994385,  450.993194580,  -20.000000000 )
create 5IRE_disulfides, elem s and bound_to elem s and /5IRE
show sphere, 5IRE_disulfides
list=[]
iterate (5IRE_disulfides), list.append((model, chain, resn, resi))
python
for i in list:
    print(i[0], i[1], i[2], i[3])
python end

quit