#150 MHz corresponds to z = 8.5
#z  freq (MHz)
#6  203
#7  178
#8  158
#9  142
#10 129
#12 109
#14 095
#16 084
#20 068
#25 055
#30 046

python mk_array_file.py -C hera350 --freq 0.203
python mk_array_file.py -C hera350 --freq 0.178
python mk_array_file.py -C hera350 --freq 0.158
python mk_array_file.py -C hera350 --freq 0.142
python mk_array_file.py -C hera350 --freq 0.129
python mk_array_file.py -C hera350 --freq 0.109
python mk_array_file.py -C hera350 --freq 0.095
python mk_array_file.py -C hera350 --freq 0.084
python mk_array_file.py -C hera350 --freq 0.068
python mk_array_file.py -C hera350 --freq 0.055
python mk_array_file.py -C hera350 --freq 0.046

python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.203GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.178GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.158GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.142GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.129GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.109GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.095GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.084GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.068GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.055GHz_arrayfile.npz -m mod
python calc_kcoverage.py hera350.drift_blmin0_blmax438_0.046GHz_arrayfile.npz -m mod
