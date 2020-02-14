from a2mdlib.utils import get_normal_modes

with open('C:/Users/Bruno/Dropbox/doctorado/aAMDlib/inp/vib/vib_modes_raw.txt') as f:
    vib = [i.strip() for i in f.readlines()]

freq, norm = get_normal_modes(
    vib, 5
)
print("HOLA")