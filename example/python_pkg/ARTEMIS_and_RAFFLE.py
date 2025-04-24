from time import sleep
import numpy
from ase.io import read, write
from ase.visualize import view
atoms = read("structures.traj", index=":")
from artemis.generator import artemis_generator
from raffle.generator import raffle_generator
art_gen = artemis_generator()
raff_gen = raffle_generator()
location, axis = art_gen.get_interface_location(atoms[0], return_fractional=True)
print("location", location)
print("axis", axis)
raff_gen.set_host(atoms[0])

missing_stoich = raff_gen.prepare_host(interface_location=location, depth=2, location_as_fractional=True)#11.97])
print("missing_stoich", missing_stoich)
host_1 = raff_gen.get_host()
view(host_1)

raff_gen.set_host(atoms[0])

host_a = raff_gen.get_host()
view(host_a)

location, axis = art_gen.get_interface_location(atoms[0], return_fractional=False)
print("location", location)
print("axis", axis)
raff_gen.set_host(atoms[0])
missing_stoich = raff_gen.prepare_host(interface_location=location, depth=3, location_as_fractional=False)#11.97])
print("missing_stoich", missing_stoich)
host_2 = raff_gen.get_host()
view(host_2)
