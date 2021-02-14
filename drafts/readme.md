NON working drafts for trying out FiPy tools

eq = (TransientTerm() == DiffusionTerm(coeff=diffCoeff) + PowerLawConvectionTerm(coeff=convCoeff) + sourceCoeff)

Copy1, Copy2 próbálkozások

Copy3:
- Dinamikus listába mentem a megoldott CellVariable objecteket, numpy arrayben is próbáltam
- ha csak egyszer szerepel az eq.solve(), akkor más eredményeket kapok, mintha többször
- egyes lépések után nem tudom kimenteni a listába a megfelelő CellVariable-eket, mert a dinamikus lista összes tagjában végül a legutolsó megoldás kerül bele
- while ciklusban sem működik
- plotolás sem megy 2D-ben (idő - CellVariable)
- legfontosabb lenne, hogy külön meglegyenek az egyes megoldások
- lehetne kimenteni fájlba, de az meghosszabbítaná a futásidőt
- lehetne 3D listát vagy arrayt létrehozni, melyben idő, sugár, és sűrűségértékek vannak de szerintem ez bonyolítana
- a CellVariable() objectnek a .value attribútumában a függvényértékek float64 arrayben vannak, de ha azt használom, akkor a mesh helykordinátái elvesznek
- DummySolver()-el dolgoztam, melynek dokumentációjában az szerepel, hogy "doesn't do anything"??!!
link: https://www.ctcms.nist.gov/fipy/fipy/generated/fipy.solvers.petsc.html#module-fipy.solvers.petsc.dummySolver
