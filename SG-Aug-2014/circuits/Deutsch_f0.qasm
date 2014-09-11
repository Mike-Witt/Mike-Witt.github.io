#
# Deutsch's algorithm, circuit for "set to zero" operator
#
# f operates on two bits, starting with bit 0
    defbox	f,2,0,'Uf_0'
    def     l01,0,'\frac{1}{\sqrt{2}}|0\rangle+\frac{1}{\sqrt{2}}|1\rangle'
    def     l02,0,'\frac{1}{\sqrt{2}}|0\rangle+\frac{1}{\sqrt{2}}|1\rangle'
    def     l11,0,'\frac{1}{\sqrt{2}}|0\rangle-\frac{1}{\sqrt{2}}|1\rangle'
    def     l12,0,'\frac{1}{\sqrt{2}}|0\rangle-\frac{1}{\sqrt{2}}|1\rangle'
    def     l03,0,'|0\rangle'
  
    qubit	0
    qubit	1

    h	0
    h	1
    l01 0
    l11 1
    f	0,1
    l02 0
    l12 1
    h	0
    l03 0
    measure 0

