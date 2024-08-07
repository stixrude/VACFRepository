        parameter (ntable=36)
        character*2 atom
        character*5 atomtab(ntable)
        real z(ntable)
	data atomtab/
     &               'H',
     &               'He',
     &               'Li',
     &               'Be',
     &               'B',
     &               'C',
     &               'N',
     &               'O',
     &               'F',
     &               'Ne',
     &               'Na',
     &               'Mg',
     &               'Al',
     &               'Si',
     &               'P',
     &               'S',
     &               'Cl',
     &               'Ar',
     &               'K',
     &               'Ca',
     &               'Sc',
     &               'Ti',
     &               'V',
     &               'Cr',
     &               'Mn',
     &               'Fe',
     &               'Co',
     &               'Ni',
     &               'Cu',
     &               'Zn',
     &               'Ga',
     &               'Ge',
     &               'As',
     &               'Se',
     &               'Br',
     &               'Kr'
     &              /
	data z/
     &               +1.,
     &               0.,
     &               +1.,
     &               +2.,
     &               +3.,
     &               +4.,
     &               +5.,
     &               -2.,
     &               -1.,
     &               0.,
     &               +1.,
     &               +2.,
     &               +3.,
     &               +4.,
     &               +5.,
     &               -2.,
     &               -1.,
     &               0.,
     &               +1.,
     &               +2.,
     &               +3.,
     &               +4.,
     &               +5.,
     &               +3.,
     &               +2.,
     &               +2.,
     &               +2.,
     &               +2.,
     &               +2.,
     &               +2.,
     &               +3.,
     &               +4.,
     &               +5.,
     &               -2.,
     &               -1.,
     &               0.
     &        /

