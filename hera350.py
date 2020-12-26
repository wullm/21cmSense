import aipy as a, numpy as n, os
import math

class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.array_params = {}
    def get_ant_params(self, ant_prms={'*':'*'}):
        prms = a.fit.AntennaArray.get_params(self, ant_prms)
        for k in ant_prms:
            top_pos = n.dot(self._eq2zen, self[int(k)].pos)
            if ant_prms[k] == '*':
                prms[k].update({'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]})
            else:
                for val in ant_prms[k]:
                    if   val == 'top_x': prms[k]['top_x'] = top_pos[0]
                    elif val == 'top_y': prms[k]['top_y'] = top_pos[1]
                    elif val == 'top_z': prms[k]['top_z'] = top_pos[2]
        return prms
    def set_ant_params(self, prms):
        changed = a.fit.AntennaArray.set_params(self, prms)
        for i, ant in enumerate(self):
            ant_changed = False
            top_pos = n.dot(self._eq2zen, ant.pos)
            try:
                top_pos[0] = prms[str(i)]['top_x']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[1] = prms[str(i)]['top_y']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[2] = prms[str(i)]['top_z']
                ant_changed = True
            except(KeyError): pass
            if ant_changed: ant.pos = n.dot(n.linalg.inv(self._eq2zen), top_pos)
            changed |= ant_changed
        return changed
    def get_arr_params(self):
        return self.array_params
    def set_arr_params(self, prms):
        for param in prms:
            self.array_params[param] = prms[param]
            if param == 'dish_size_in_lambda':
                FWHM = 2.35*(0.45/prms[param]) #radians
                self.array_params['obs_duration'] = 61.*FWHM / (15.*a.const.deg)# minutes it takes the sky to drift through beam FWHM
            if param == 'antpos':
                bl_lens = n.sum(n.array(prms[param])**2,axis=1)**.5
        return self.array_params

#===========================ARRAY SPECIFIC PARAMETERS==========================

#Set antenna positions here; for regular arrays like Hera we can use an algorithm; otherwise antpos should just be a list of [x,y,z] coords in light-nanoseconds
##############################
##############################
SScale = 1.0
Dant = 1400.0*SScale  # cm
Spacing = 60.0*SScale # cm
nside = 11.           # hex number
freq = 150.0          # MHz
##############################
##############################

L = (Dant+Spacing) / a.const.len_ns
dL = (Dant+Spacing)*math.cos(30.0*math.pi/180.0) / a.const.len_ns #close packed hex
ant_size_wavel = (Dant/100.0) / (300./freq)

# antpos = []
# cen_y, cen_z = 0, 0
# for row in n.arange(nside):
#     for cen_x in n.arange((2*nside-1)-row):
#         dx = row/2
#         antpos.append(((cen_x + dx)*L, row*dL, cen_z))
#         if row != 0:
#             antpos.append(((cen_x + dx)*L, -row*dL, cen_z))

antpos = [ (-8.029999999999999716e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(-6.570000000000000284e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(-5.109999999999999432e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(-3.650000000000000000e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(-2.189999999999999858e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(-7.299999999999999822e+00, 1.306543659176122958e+02, 0.000000000000000000e+00),
(7.299999999999999822e+00, 1.306543659176122958e+02, 0.000000000000000000e+00),
(2.189999999999999858e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(3.650000000000000000e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(5.110000000000000142e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(6.570000000000000284e+01, 1.306543659176122958e+02, 0.000000000000000000e+00),
(-8.759999999999999432e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(-7.300000000000000000e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(-5.839999999999999858e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(-4.379999999999999716e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(-2.919999999999999929e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(-1.459999999999999964e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(0.000000000000000000e+00, 1.180103950223595035e+02, 0.000000000000000000e+00),
(1.459999999999999787e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(2.919999999999999929e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(4.380000000000000426e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(5.840000000000000568e+01, 1.180103950223595035e+02, 0.000000000000000000e+00),
(8.029999999999999716e+01, 1.222250519874437629e+02, 0.000000000000000000e+00),
(-9.489999999999999147e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(-8.029999999999999716e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(-6.570000000000000284e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(-5.109999999999999432e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(-3.650000000000000000e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(-2.189999999999999858e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(-7.299999999999999822e+00, 1.053664241271066970e+02, 0.000000000000000000e+00),
(7.299999999999999822e+00, 1.053664241271066970e+02, 0.000000000000000000e+00),
(2.189999999999999858e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(3.650000000000000000e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(5.110000000000000142e+01, 1.053664241271066970e+02, 0.000000000000000000e+00),
(7.300000000000000000e+01, 1.095810810921909564e+02, 0.000000000000000000e+00),
(8.759999999999999432e+01, 1.095810810921909564e+02, 0.000000000000000000e+00),
(-1.021999999999999886e+02, 9.272245323185390475e+01, 0.000000000000000000e+00),
(-8.759999999999999432e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(-7.300000000000000000e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(-5.839999999999999858e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(-4.379999999999999716e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(-2.919999999999999929e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(-1.459999999999999964e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(0.000000000000000000e+00, 9.272245323185390475e+01, 0.000000000000000000e+00),
(1.459999999999999787e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(2.919999999999999929e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(4.380000000000000426e+01, 9.272245323185390475e+01, 0.000000000000000000e+00),
(6.570000000000000284e+01, 9.693711019693816411e+01, 0.000000000000000000e+00),
(8.029999999999999716e+01, 9.693711019693816411e+01, 0.000000000000000000e+00),
(9.489999999999999147e+01, 9.693711019693816411e+01, 0.000000000000000000e+00),
(-1.095000000000000000e+02, 8.007848233660108406e+01, 0.000000000000000000e+00),
(-9.489999999999999147e+01, 8.007848233660108406e+01, 0.000000000000000000e+00),
(-8.029999999999999716e+01, 8.007848233660108406e+01, 0.000000000000000000e+00),
(-6.570000000000000284e+01, 8.007848233660108406e+01, 0.000000000000000000e+00),
(-5.109999999999999432e+01, 8.007848233660108406e+01, 0.000000000000000000e+00),
(-3.650000000000000000e+01, 8.007848233660108406e+01, 0.000000000000000000e+00),
(-2.189999999999999858e+01, 8.007848233660108406e+01, 0.000000000000000000e+00),
(-7.299999999999999822e+00, 8.007848233660108406e+01, 0.000000000000000000e+00),
(7.299999999999999822e+00, 8.007848233660108406e+01, 0.000000000000000000e+00),
(2.189999999999999858e+01, 8.007848233660108406e+01, 0.000000000000000000e+00),
(3.650000000000000000e+01, 8.007848233660108406e+01, 0.000000000000000000e+00),
(5.839999999999999858e+01, 8.429313930168534341e+01, 0.000000000000000000e+00),
(7.300000000000000000e+01, 8.429313930168534341e+01, 0.000000000000000000e+00),
(8.759999999999999432e+01, 8.429313930168534341e+01, 0.000000000000000000e+00),
(1.022000000000000028e+02, 8.429313930168534341e+01, 0.000000000000000000e+00),
(-1.167999999999999972e+02, 6.743451144134829178e+01, 0.000000000000000000e+00),
(-1.021999999999999886e+02, 6.743451144134829178e+01, 0.000000000000000000e+00),
(-8.759999999999999432e+01, 6.743451144134829178e+01, 0.000000000000000000e+00),
(-7.300000000000000000e+01, 6.743451144134829178e+01, 0.000000000000000000e+00),
(-5.839999999999999858e+01, 6.743451144134829178e+01, 0.000000000000000000e+00),
(-4.379999999999999716e+01, 6.743451144134829178e+01, 0.000000000000000000e+00),
(-2.919999999999999929e+01, 6.743451144134829178e+01, 0.000000000000000000e+00),
(-1.459999999999999964e+01, 6.743451144134829178e+01, 0.000000000000000000e+00),
(0.000000000000000000e+00, 6.743451144134829178e+01, 0.000000000000000000e+00),
(1.459999999999999787e+01, 6.743451144134829178e+01, 0.000000000000000000e+00),
(2.919999999999999929e+01, 6.743451144134829178e+01, 0.000000000000000000e+00),
(5.110000000000000142e+01, 7.164916840643255114e+01, 0.000000000000000000e+00),
(6.570000000000000284e+01, 7.164916840643255114e+01, 0.000000000000000000e+00),
(8.029999999999999716e+01, 7.164916840643255114e+01, 0.000000000000000000e+00),
(9.489999999999999147e+01, 7.164916840643255114e+01, 0.000000000000000000e+00),
(1.095000000000000000e+02, 7.164916840643255114e+01, 0.000000000000000000e+00),
(-1.240999999999999943e+02, 5.479054054609548530e+01, 0.000000000000000000e+00),
(-1.095000000000000000e+02, 5.479054054609548530e+01, 0.000000000000000000e+00),
(-9.489999999999999147e+01, 5.479054054609548530e+01, 0.000000000000000000e+00),
(-8.029999999999999716e+01, 5.479054054609548530e+01, 0.000000000000000000e+00),
(-6.570000000000000284e+01, 5.479054054609548530e+01, 0.000000000000000000e+00),
(-5.109999999999999432e+01, 5.479054054609548530e+01, 0.000000000000000000e+00),
(-3.650000000000000000e+01, 5.479054054609548530e+01, 0.000000000000000000e+00),
(-2.189999999999999858e+01, 5.479054054609548530e+01, 0.000000000000000000e+00),
(-7.299999999999999822e+00, 5.479054054609548530e+01, 0.000000000000000000e+00),
(7.299999999999999822e+00, 5.479054054609548530e+01, 0.000000000000000000e+00),
(2.189999999999999858e+01, 5.479054054609548530e+01, 0.000000000000000000e+00),
(4.379999999999999716e+01, 5.900519751117974465e+01, 0.000000000000000000e+00),
(5.839999999999999858e+01, 5.900519751117974465e+01, 0.000000000000000000e+00),
(7.300000000000000000e+01, 5.900519751117974465e+01, 0.000000000000000000e+00),
(8.759999999999999432e+01, 5.900519751117974465e+01, 0.000000000000000000e+00),
(1.022000000000000028e+02, 5.900519751117974465e+01, 0.000000000000000000e+00),
(1.167999999999999972e+02, 5.900519751117974465e+01, 0.000000000000000000e+00),
(-1.314000000000000057e+02, 4.214656965084267881e+01, 0.000000000000000000e+00),
(-1.167999999999999972e+02, 4.214656965084267881e+01, 0.000000000000000000e+00),
(-1.021999999999999886e+02, 4.214656965084267881e+01, 0.000000000000000000e+00),
(-8.759999999999999432e+01, 4.214656965084267881e+01, 0.000000000000000000e+00),
(-7.300000000000000000e+01, 4.214656965084267881e+01, 0.000000000000000000e+00),
(-5.839999999999999858e+01, 4.214656965084267881e+01, 0.000000000000000000e+00),
(-4.379999999999999716e+01, 4.214656965084267881e+01, 0.000000000000000000e+00),
(-2.919999999999999929e+01, 4.214656965084267881e+01, 0.000000000000000000e+00),
(-1.459999999999999964e+01, 4.214656965084267881e+01, 0.000000000000000000e+00),
(0.000000000000000000e+00, 4.214656965084267881e+01, 0.000000000000000000e+00),
(1.459999999999999787e+01, 4.214656965084267881e+01, 0.000000000000000000e+00),
(3.650000000000000000e+01, 4.636122661592693817e+01, 0.000000000000000000e+00),
(5.110000000000000142e+01, 4.636122661592693817e+01, 0.000000000000000000e+00),
(6.570000000000000284e+01, 4.636122661592693817e+01, 0.000000000000000000e+00),
(8.029999999999999716e+01, 4.636122661592693817e+01, 0.000000000000000000e+00),
(9.489999999999999147e+01, 4.636122661592693817e+01, 0.000000000000000000e+00),
(1.095000000000000000e+02, 4.636122661592693817e+01, 0.000000000000000000e+00),
(1.240999999999999943e+02, 4.636122661592693817e+01, 0.000000000000000000e+00),
(-1.387000000000000171e+02, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-1.240999999999999943e+02, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-1.095000000000000000e+02, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-9.489999999999999147e+01, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-8.029999999999999716e+01, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-6.570000000000000284e+01, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-5.109999999999999432e+01, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-3.650000000000000000e+01, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-2.189999999999999858e+01, 2.950259875558987233e+01, 0.000000000000000000e+00),
(-7.299999999999999822e+00, 2.950259875558987233e+01, 0.000000000000000000e+00),
(7.299999999999999822e+00, 2.950259875558987233e+01, 0.000000000000000000e+00),
(2.919999999999999929e+01, 3.371725572067413879e+01, 0.000000000000000000e+00),
(4.379999999999999716e+01, 3.371725572067413879e+01, 0.000000000000000000e+00),
(5.839999999999999858e+01, 3.371725572067413879e+01, 0.000000000000000000e+00),
(7.300000000000000000e+01, 3.371725572067413879e+01, 0.000000000000000000e+00),
(8.759999999999999432e+01, 3.371725572067413879e+01, 0.000000000000000000e+00),
(1.022000000000000028e+02, 3.371725572067413879e+01, 0.000000000000000000e+00),
(1.167999999999999972e+02, 3.371725572067413879e+01, 0.000000000000000000e+00),
(1.314000000000000057e+02, 3.371725572067413879e+01, 0.000000000000000000e+00),
(-1.460000000000000000e+02, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-1.314000000000000057e+02, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-1.167999999999999972e+02, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-1.021999999999999886e+02, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-8.759999999999999432e+01, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-7.300000000000000000e+01, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-5.839999999999999858e+01, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-4.379999999999999716e+01, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-2.919999999999999929e+01, 1.685862786033707295e+01, 0.000000000000000000e+00),
(-1.459999999999999964e+01, 1.685862786033707295e+01, 0.000000000000000000e+00),
(0.000000000000000000e+00, 1.685862786033707295e+01, 0.000000000000000000e+00),
(2.189999999999999858e+01, 2.107328482542133941e+01, 0.000000000000000000e+00),
(3.650000000000000000e+01, 2.107328482542133941e+01, 0.000000000000000000e+00),
(5.110000000000000142e+01, 2.107328482542133941e+01, 0.000000000000000000e+00),
(6.570000000000000284e+01, 2.107328482542133941e+01, 0.000000000000000000e+00),
(8.029999999999999716e+01, 2.107328482542133941e+01, 0.000000000000000000e+00),
(9.489999999999999147e+01, 2.107328482542133941e+01, 0.000000000000000000e+00),
(1.095000000000000000e+02, 2.107328482542133941e+01, 0.000000000000000000e+00),
(1.240999999999999943e+02, 2.107328482542133941e+01, 0.000000000000000000e+00),
(1.386999999999999886e+02, 2.107328482542133941e+01, 0.000000000000000000e+00),
(-1.460000000000000000e+02, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-1.314000000000000057e+02, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-1.167999999999999972e+02, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-1.022000000000000028e+02, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-8.759999999999999432e+01, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-7.300000000000000000e+01, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-5.839999999999999858e+01, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-4.379999999999999716e+01, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-2.919999999999999929e+01, 0.000000000000000000e+00, 0.000000000000000000e+00),
(-1.459999999999999964e+01, 0.000000000000000000e+00, 0.000000000000000000e+00),
(0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00),
(1.459999999999999964e+01, 8.429313930168534696e+00, 0.000000000000000000e+00),
(2.919999999999999929e+01, 8.429313930168534696e+00, 0.000000000000000000e+00),
(4.379999999999999716e+01, 8.429313930168534696e+00, 0.000000000000000000e+00),
(5.839999999999999858e+01, 8.429313930168534696e+00, 0.000000000000000000e+00),
(7.300000000000000000e+01, 8.429313930168534696e+00, 0.000000000000000000e+00),
(8.759999999999999432e+01, 8.429313930168534696e+00, 0.000000000000000000e+00),
(1.022000000000000028e+02, 8.429313930168534696e+00, 0.000000000000000000e+00),
(1.167999999999999972e+02, 8.429313930168534696e+00, 0.000000000000000000e+00),
(1.314000000000000057e+02, 8.429313930168534696e+00, 0.000000000000000000e+00),
(1.460000000000000000e+02, 8.429313930168534696e+00, 0.000000000000000000e+00),
(-1.386999999999999886e+02, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-1.240999999999999943e+02, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-1.095000000000000000e+02, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-9.489999999999999147e+01, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-8.029999999999999716e+01, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-6.570000000000000284e+01, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-5.110000000000000142e+01, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-3.650000000000000000e+01, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-2.189999999999999858e+01, -1.264397089525280293e+01, 0.000000000000000000e+00),
(-7.299999999999999822e+00, -1.264397089525280293e+01, 0.000000000000000000e+00),
(7.299999999999999822e+00, -1.264397089525280293e+01, 0.000000000000000000e+00),
(2.189999999999999858e+01, -4.214656965084268236e+00, 0.000000000000000000e+00),
(3.650000000000000000e+01, -4.214656965084268236e+00, 0.000000000000000000e+00),
(5.110000000000000142e+01, -4.214656965084268236e+00, 0.000000000000000000e+00),
(6.570000000000000284e+01, -4.214656965084268236e+00, 0.000000000000000000e+00),
(8.029999999999999716e+01, -4.214656965084268236e+00, 0.000000000000000000e+00),
(9.489999999999999147e+01, -4.214656965084268236e+00, 0.000000000000000000e+00),
(1.095000000000000000e+02, -4.214656965084268236e+00, 0.000000000000000000e+00),
(1.240999999999999943e+02, -4.214656965084268236e+00, 0.000000000000000000e+00),
(1.386999999999999886e+02, -4.214656965084268236e+00, 0.000000000000000000e+00),
(-1.314000000000000057e+02, -2.528794179050560587e+01, 0.000000000000000000e+00),
(-1.167999999999999972e+02, -2.528794179050560587e+01, 0.000000000000000000e+00),
(-1.022000000000000028e+02, -2.528794179050560587e+01, 0.000000000000000000e+00),
(-8.759999999999999432e+01, -2.528794179050560587e+01, 0.000000000000000000e+00),
(-7.300000000000000000e+01, -2.528794179050560587e+01, 0.000000000000000000e+00),
(-5.839999999999999858e+01, -2.528794179050560587e+01, 0.000000000000000000e+00),
(-4.379999999999999716e+01, -2.528794179050560587e+01, 0.000000000000000000e+00),
(-2.919999999999999929e+01, -2.528794179050560587e+01, 0.000000000000000000e+00),
(-1.459999999999999964e+01, -2.528794179050560587e+01, 0.000000000000000000e+00),
(0.000000000000000000e+00, -2.528794179050560587e+01, 0.000000000000000000e+00),
(1.459999999999999964e+01, -2.528794179050560587e+01, 0.000000000000000000e+00),
(2.919999999999999929e+01, -1.685862786033707295e+01, 0.000000000000000000e+00),
(4.379999999999999716e+01, -1.685862786033707295e+01, 0.000000000000000000e+00),
(5.839999999999999858e+01, -1.685862786033707295e+01, 0.000000000000000000e+00),
(7.300000000000000000e+01, -1.685862786033707295e+01, 0.000000000000000000e+00),
(8.759999999999999432e+01, -1.685862786033707295e+01, 0.000000000000000000e+00),
(1.022000000000000028e+02, -1.685862786033707295e+01, 0.000000000000000000e+00),
(1.167999999999999972e+02, -1.685862786033707295e+01, 0.000000000000000000e+00),
(1.314000000000000057e+02, -1.685862786033707295e+01, 0.000000000000000000e+00),
(-1.240999999999999943e+02, -3.793191268575840525e+01, 0.000000000000000000e+00),
(-1.095000000000000000e+02, -3.793191268575840525e+01, 0.000000000000000000e+00),
(-9.489999999999999147e+01, -3.793191268575840525e+01, 0.000000000000000000e+00),
(-8.029999999999999716e+01, -3.793191268575840525e+01, 0.000000000000000000e+00),
(-6.570000000000000284e+01, -3.793191268575840525e+01, 0.000000000000000000e+00),
(-5.110000000000000142e+01, -3.793191268575840525e+01, 0.000000000000000000e+00),
(-3.650000000000000000e+01, -3.793191268575840525e+01, 0.000000000000000000e+00),
(-2.189999999999999858e+01, -3.793191268575840525e+01, 0.000000000000000000e+00),
(-7.299999999999999822e+00, -3.793191268575840525e+01, 0.000000000000000000e+00),
(7.299999999999999822e+00, -3.793191268575840525e+01, 0.000000000000000000e+00),
(2.189999999999999858e+01, -3.793191268575840525e+01, 0.000000000000000000e+00),
(3.650000000000000000e+01, -2.950259875558987233e+01, 0.000000000000000000e+00),
(5.110000000000000142e+01, -2.950259875558987233e+01, 0.000000000000000000e+00),
(6.570000000000000284e+01, -2.950259875558987233e+01, 0.000000000000000000e+00),
(8.029999999999999716e+01, -2.950259875558987233e+01, 0.000000000000000000e+00),
(9.489999999999999147e+01, -2.950259875558987233e+01, 0.000000000000000000e+00),
(1.095000000000000000e+02, -2.950259875558987233e+01, 0.000000000000000000e+00),
(1.240999999999999943e+02, -2.950259875558987233e+01, 0.000000000000000000e+00),
(-1.167999999999999972e+02, -5.057588358101121173e+01, 0.000000000000000000e+00),
(-1.022000000000000028e+02, -5.057588358101121173e+01, 0.000000000000000000e+00),
(-8.759999999999999432e+01, -5.057588358101121173e+01, 0.000000000000000000e+00),
(-7.300000000000000000e+01, -5.057588358101121173e+01, 0.000000000000000000e+00),
(-5.839999999999999858e+01, -5.057588358101121173e+01, 0.000000000000000000e+00),
(-4.379999999999999716e+01, -5.057588358101121173e+01, 0.000000000000000000e+00),
(-2.919999999999999929e+01, -5.057588358101121173e+01, 0.000000000000000000e+00),
(-1.459999999999999964e+01, -5.057588358101121173e+01, 0.000000000000000000e+00),
(0.000000000000000000e+00, -5.057588358101121173e+01, 0.000000000000000000e+00),
(1.459999999999999964e+01, -5.057588358101121173e+01, 0.000000000000000000e+00),
(2.919999999999999929e+01, -5.057588358101121173e+01, 0.000000000000000000e+00),
(4.379999999999999716e+01, -4.214656965084267881e+01, 0.000000000000000000e+00),
(5.839999999999999858e+01, -4.214656965084267881e+01, 0.000000000000000000e+00),
(7.300000000000000000e+01, -4.214656965084267881e+01, 0.000000000000000000e+00),
(8.759999999999999432e+01, -4.214656965084267881e+01, 0.000000000000000000e+00),
(1.022000000000000028e+02, -4.214656965084267881e+01, 0.000000000000000000e+00),
(1.167999999999999972e+02, -4.214656965084267881e+01, 0.000000000000000000e+00),
(-1.095000000000000000e+02, -6.321985447626401822e+01, 0.000000000000000000e+00),
(-9.489999999999999147e+01, -6.321985447626401822e+01, 0.000000000000000000e+00),
(-8.029999999999999716e+01, -6.321985447626401822e+01, 0.000000000000000000e+00),
(-6.570000000000000284e+01, -6.321985447626401822e+01, 0.000000000000000000e+00),
(-5.110000000000000142e+01, -6.321985447626401822e+01, 0.000000000000000000e+00),
(-3.650000000000000000e+01, -6.321985447626401822e+01, 0.000000000000000000e+00),
(-2.189999999999999858e+01, -6.321985447626401822e+01, 0.000000000000000000e+00),
(-7.299999999999999822e+00, -6.321985447626401822e+01, 0.000000000000000000e+00),
(7.299999999999999822e+00, -6.321985447626401822e+01, 0.000000000000000000e+00),
(2.189999999999999858e+01, -6.321985447626401822e+01, 0.000000000000000000e+00),
(3.650000000000000000e+01, -6.321985447626401822e+01, 0.000000000000000000e+00),
(5.110000000000000142e+01, -5.479054054609548530e+01, 0.000000000000000000e+00),
(6.570000000000000284e+01, -5.479054054609548530e+01, 0.000000000000000000e+00),
(8.029999999999999716e+01, -5.479054054609548530e+01, 0.000000000000000000e+00),
(9.489999999999999147e+01, -5.479054054609548530e+01, 0.000000000000000000e+00),
(1.095000000000000000e+02, -5.479054054609548530e+01, 0.000000000000000000e+00),
(-1.022000000000000028e+02, -7.586382537151681049e+01, 0.000000000000000000e+00),
(-8.759999999999999432e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(-7.300000000000000000e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(-5.839999999999999858e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(-4.379999999999999716e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(-2.919999999999999929e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(-1.459999999999999964e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(0.000000000000000000e+00, -7.586382537151681049e+01, 0.000000000000000000e+00),
(1.459999999999999964e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(2.919999999999999929e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(4.379999999999999716e+01, -7.586382537151681049e+01, 0.000000000000000000e+00),
(5.839999999999999858e+01, -6.743451144134827757e+01, 0.000000000000000000e+00),
(7.300000000000000000e+01, -6.743451144134827757e+01, 0.000000000000000000e+00),
(8.759999999999999432e+01, -6.743451144134827757e+01, 0.000000000000000000e+00),
(1.022000000000000028e+02, -6.743451144134827757e+01, 0.000000000000000000e+00),
(-9.489999999999999147e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(-8.029999999999999716e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(-6.570000000000000284e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(-5.110000000000000142e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(-3.650000000000000000e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(-2.189999999999999858e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(-7.299999999999999822e+00, -8.850779626676963119e+01, 0.000000000000000000e+00),
(7.299999999999999822e+00, -8.850779626676963119e+01, 0.000000000000000000e+00),
(2.189999999999999858e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(3.650000000000000000e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(5.110000000000000142e+01, -8.850779626676963119e+01, 0.000000000000000000e+00),
(6.570000000000000284e+01, -8.007848233660109827e+01, 0.000000000000000000e+00),
(8.029999999999999716e+01, -8.007848233660109827e+01, 0.000000000000000000e+00),
(9.489999999999999147e+01, -8.007848233660109827e+01, 0.000000000000000000e+00),
(-8.759999999999999432e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(-7.300000000000000000e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(-5.839999999999999858e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(-4.379999999999999716e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(-2.919999999999999929e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(-1.459999999999999964e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(0.000000000000000000e+00, -1.011517671620224235e+02, 0.000000000000000000e+00),
(1.459999999999999964e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(2.919999999999999929e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(4.379999999999999716e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(5.839999999999999858e+01, -1.011517671620224235e+02, 0.000000000000000000e+00),
(7.300000000000000000e+01, -9.272245323185389054e+01, 0.000000000000000000e+00),
(8.759999999999999432e+01, -9.272245323185389054e+01, 0.000000000000000000e+00),
(-8.029999999999999716e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(-6.570000000000000284e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(-5.110000000000000142e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(-3.650000000000000000e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(-2.189999999999999858e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(-7.299999999999999822e+00, -1.137957380572752299e+02, 0.000000000000000000e+00),
(7.299999999999999822e+00, -1.137957380572752299e+02, 0.000000000000000000e+00),
(2.189999999999999858e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(3.650000000000000000e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(5.110000000000000142e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(6.570000000000000284e+01, -1.137957380572752299e+02, 0.000000000000000000e+00),
(8.029999999999999716e+01, -1.053664241271066970e+02, 0.000000000000000000e+00),
(-2.190000000000000000e+02, 3.456018711369099492e+02, 0.000000000000000000e+00),
(-7.300000000000000000e+01, 3.456018711369099492e+02, 0.000000000000000000e+00),
(7.300000000000000000e+01, 3.456018711369099492e+02, 0.000000000000000000e+00),
(2.190000000000000000e+02, 3.456018711369099492e+02, 0.000000000000000000e+00),
(-2.920000000000000000e+02, 2.275914761145504599e+02, 0.000000000000000000e+00),
(-1.460000000000000000e+02, 2.191621621843819412e+02, 0.000000000000000000e+00),
(0.000000000000000000e+00, 2.191621621843819412e+02, 0.000000000000000000e+00),
(1.460000000000000000e+02, 2.191621621843819412e+02, 0.000000000000000000e+00),
(2.920000000000000000e+02, 2.191621621843819412e+02, 0.000000000000000000e+00),
(-3.650000000000000000e+02, 1.011517671620224377e+02, 0.000000000000000000e+00),
(-2.190000000000000000e+02, 1.011517671620224377e+02, 0.000000000000000000e+00),
(2.190000000000000000e+02, 9.272245323185390475e+01, 0.000000000000000000e+00),
(3.650000000000000000e+02, 9.272245323185390475e+01, 0.000000000000000000e+00),
(-4.380000000000000000e+02, -2.528794179050560231e+01, 0.000000000000000000e+00),
(-2.920000000000000000e+02, -2.528794179050560231e+01, 0.000000000000000000e+00),
(2.920000000000000000e+02, -1.685862786033706939e+01, 0.000000000000000000e+00),
(4.380000000000000000e+02, -1.685862786033706939e+01, 0.000000000000000000e+00),
(-3.650000000000000000e+02, -1.517276507430336494e+02, 0.000000000000000000e+00),
(-2.190000000000000000e+02, -1.517276507430336494e+02, 0.000000000000000000e+00),
(2.190000000000000000e+02, -1.432983368128651023e+02, 0.000000000000000000e+00),
(3.650000000000000000e+02, -1.432983368128651023e+02, 0.000000000000000000e+00),
(-2.920000000000000000e+02, -2.781673596955616858e+02, 0.000000000000000000e+00),
(-1.460000000000000000e+02, -2.781673596955616858e+02, 0.000000000000000000e+00),
(0.000000000000000000e+00, -2.697380457653931671e+02, 0.000000000000000000e+00),
(1.460000000000000000e+02, -2.697380457653931671e+02, 0.000000000000000000e+00),
(2.920000000000000000e+02, -2.697380457653931671e+02, 0.000000000000000000e+00),
(-2.190000000000000000e+02, -4.046070686480896939e+02, 0.000000000000000000e+00),
(-7.300000000000000000e+01, -3.961777547179211751e+02, 0.000000000000000000e+00),
(7.300000000000000000e+01, -3.961777547179211751e+02, 0.000000000000000000e+00),
(2.190000000000000000e+02, -3.961777547179211751e+02, 0.000000000000000000e+00)]

#Convert from metres to light-nanoseconds
for i,pos in enumerate(antpos):
    antpos[i] = n.array(pos) * 100 / a.const.len_ns

#Set other array parameters here
prms = {
    'name': os.path.basename(__file__)[:-3], #remove .py from filename
    'loc': ('-30:43:17.5', '21:25:41.9'), # Karoo  #'loc': ('38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
    'antpos': antpos,
    'beam': a.fit.Beam2DGaussian,
    'dish_size_in_lambda': ant_size_wavel, #in units of wavelengths at 150 MHz = 2 meters; this will also define the observation duration
    'Trx': 1e5 #receiver temp in mK
}

#=======================END ARRAY SPECIFIC PARAMETERS==========================

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        beam = prms['beam'](freqs, xwidth=(0.45/prms['dish_size_in_lambda']), ywidth=(0.45/prms['dish_size_in_lambda'])) #as it stands, the size of the beam as defined here is not actually used anywhere in this package, but is a necessary parameter for the aipy Beam2DGaussian object
        antennas.append(a.fit.Antenna(0, 0, 0, beam))
    aa = AntennaArray(prms['loc'], antennas)
    p = {}
    for i in range(nants):
        top_pos = prms['antpos'][i]
        p[str(i)] = {'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]}
    aa.set_ant_params(p)
    aa.set_arr_params(prms)
    return aa

def get_catalog(*args, **kwargs): return a.src.get_catalog(*args, **kwargs)