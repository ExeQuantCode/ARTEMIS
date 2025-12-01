from ase import Atoms
from ase.build import bulk
from ase.io import write
from artemis.generator import artemis_generator

generator = artemis_generator()

Si = bulk('Si', 'diamond', a=5.43, cubic=True)
Ge = bulk('Ge', 'diamond', a=5.66, cubic=True)

generator.set_materials(Si, Ge)

generator.set_surface_properties(
    miller_lw = [ 1, 1, 0 ],
    miller_up = [ 1, 1, 0 ],
)

generator.set_shift_method(num_shifts = 1)
generator.set_match_method(max_num_matches = 1)
structures = generator.generate(verbose=1)#calc=calc)


write('structures.traj', structures)

output = generator.get_all_structures_data()
print(output)

output = generator.get_structure_data(0)
print(output)