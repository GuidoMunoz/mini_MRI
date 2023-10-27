T = 56000e-6 # Tiempo del pulso en total (s)
N = 32000 #Número de muestras
Ns = T/N
#Pulso de 90 grados dura 80 us, y de 180 grados dura 160 us
#Pulso de desfase
init = 200e-6 #Pendiente de bajada
start_neg = 450e-6 #Planicie
neg = 1450e-6 #Pendiente de subida
stop_neg = 1700e-6
#Pulso de refase
init_refase = 2570e-6 #--0.12
start_pos = 3070e-6 #Pendiente de subida
pos = 5070e-6 #Pos
end_pos = 5570e-6 #Pendiente de bajada
#Spoiler
start_spoiler = 6000e-6
spoiler = 12000e-6
stop_spoiler = 50000e-6
end_spoiler = 56000e-6
v_high = 1
v_low = 1
A = 1
B = v_low/v_high 
C = 1
#######
f = open("grads0709_spoiler_real.txt", "w")
for i in range(1,N):
    if i*Ns < init:
        f.write(str(0)+"\n") 
    elif i*Ns < start_neg:
        m = (A/(start_neg-init))
        int_zero = m*init
        val =  m*Ns*i - int_zero
        f.write(str(val)+"\n") 
    elif i*Ns < neg:
        f.write(str(A)+"\n")
    elif i*Ns < stop_neg:
        m = A/(stop_neg - neg)
        int_zero = -m*stop_neg
        val = m*Ns*i + int_zero
        f.write(str(-val)+"\n")
    elif i*Ns < init_refase:
        f.write(str(0)+"\n") 
    elif i*Ns < start_pos:
        m = B/(start_pos - init_refase)
        int_zero = m*init_refase
        val =  m*Ns*i - int_zero
        f.write(str(val)+"\n")
    elif i*Ns < pos:
        f.write(str(B)+"\n")
    elif i*Ns < end_pos:
        m = B/(end_pos- pos)
        int_zero = m*end_pos
        val = m*Ns*i - int_zero
        f.write(str(-val)+"\n")
    elif i*Ns < start_spoiler:
        f.write(str(0)+"\n")
    elif i*Ns < spoiler:
        m = C/(spoiler - start_spoiler)
        int_zero = m*start_spoiler
        val =  m*Ns*i - int_zero
        f.write(str(val)+"\n")
    elif i*Ns < stop_spoiler:
        f.write(str(C)+"\n")
    elif i*Ns < end_spoiler:
        m = C/(end_spoiler- stop_spoiler)
        int_zero = m*end_spoiler
        val = m*Ns*i - int_zero
        f.write(str(-val)+"\n")
    else:
        f.write(str(0)+"\n")
print("archivo_escrito_correctamente")