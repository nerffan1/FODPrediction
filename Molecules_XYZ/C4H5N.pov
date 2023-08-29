#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {orthographic
  right -7.13*x up 7.13*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}


light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}
// no fog
#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.070;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     translate LOC}
#end

// no cell vertices
atom(< -1.32,  -2.81,  -2.20>, 0.42, rgb <0.18, 0.31, 0.97>, 0.0, ase3) // #0
atom(< -0.45,  -2.01,  -2.13>, 0.44, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #1
atom(<  0.60,  -1.05,  -2.06>, 0.44, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #2
atom(< -0.84,   0.78,  -1.05>, 0.44, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #3
atom(<  0.44,   0.20,  -1.57>, 0.44, rgb <0.56, 0.56, 0.56>, 0.0, ase3) // #4
atom(<  1.32,   0.85,  -1.56>, 0.18, rgb <1.00, 1.00, 1.00>, 0.0, ase3) // #5
atom(< -1.67,   0.06,  -1.10>, 0.18, rgb <1.00, 1.00, 1.00>, 0.0, ase3) // #6
atom(< -1.11,   1.68,  -1.62>, 0.18, rgb <1.00, 1.00, 1.00>, 0.0, ase3) // #7
atom(< -0.71,   1.11,   0.00>, 0.18, rgb <1.00, 1.00, 1.00>, 0.0, ase3) // #8
atom(<  1.58,  -1.39,  -2.42>, 0.18, rgb <1.00, 1.00, 1.00>, 0.0, ase3) // #9
cylinder {< -1.32,  -2.81,  -2.20>, < -0.89,  -2.41,  -2.17>, Rbond texture{pigment {color rgb <0.18, 0.31, 0.97> transmit 0.0} finish{ase3}}}
cylinder {< -0.45,  -2.01,  -2.13>, < -0.89,  -2.41,  -2.17>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {< -0.45,  -2.01,  -2.13>, <  0.08,  -1.53,  -2.09>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {<  0.60,  -1.05,  -2.06>, <  0.08,  -1.53,  -2.09>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {<  0.60,  -1.05,  -2.06>, <  0.52,  -0.43,  -1.81>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {<  0.44,   0.20,  -1.57>, <  0.52,  -0.43,  -1.81>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {<  0.60,  -1.05,  -2.06>, <  1.09,  -1.22,  -2.24>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {<  1.58,  -1.39,  -2.42>, <  1.09,  -1.22,  -2.24>, Rbond texture{pigment {color rgb <1.00, 1.00, 1.00> transmit 0.0} finish{ase3}}}
cylinder {< -0.84,   0.78,  -1.05>, < -0.20,   0.49,  -1.31>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {<  0.44,   0.20,  -1.57>, < -0.20,   0.49,  -1.31>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {< -0.84,   0.78,  -1.05>, < -1.25,   0.42,  -1.07>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {< -1.67,   0.06,  -1.10>, < -1.25,   0.42,  -1.07>, Rbond texture{pigment {color rgb <1.00, 1.00, 1.00> transmit 0.0} finish{ase3}}}
cylinder {< -0.84,   0.78,  -1.05>, < -0.98,   1.23,  -1.34>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {< -1.11,   1.68,  -1.62>, < -0.98,   1.23,  -1.34>, Rbond texture{pigment {color rgb <1.00, 1.00, 1.00> transmit 0.0} finish{ase3}}}
cylinder {< -0.84,   0.78,  -1.05>, < -0.77,   0.94,  -0.52>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {< -0.71,   1.11,   0.00>, < -0.77,   0.94,  -0.52>, Rbond texture{pigment {color rgb <1.00, 1.00, 1.00> transmit 0.0} finish{ase3}}}
cylinder {<  0.44,   0.20,  -1.57>, <  0.88,   0.52,  -1.56>, Rbond texture{pigment {color rgb <0.56, 0.56, 0.56> transmit 0.0} finish{ase3}}}
cylinder {<  1.32,   0.85,  -1.56>, <  0.88,   0.52,  -1.56>, Rbond texture{pigment {color rgb <1.00, 1.00, 1.00> transmit 0.0} finish{ase3}}}
// no constraints
