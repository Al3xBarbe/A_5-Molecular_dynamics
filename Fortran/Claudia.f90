!========================== generar_normal(s, a0, c0, m0, no1, no2, In2) ===========================================
!========================== generar_normal_txt_congr(pasos, semilla, a, c, m, nombre_archivo) ===========================================
!========================== estudio_rand_f90(iteraciones, nombre_archivo) ===========================================
!========================== generar_normal_txt_f90(pasos, nombre_archivo) ===========================================
!========================== generar_normal_f90(q1,q2) ===========================================
module generador_normal
    contains 
    subroutine generar_normal(s, a0, c0, m0, no1, no2, In2)

        implicit none

        integer(kind=8), intent(in) :: s, a0, c0, m0
        integer(kind=8), intent(out) :: In2
        
        real, intent(out) :: no1, no2

        integer (kind=8) :: I1, I2
        real :: x, y, pi

        pi = 4.0 * atan(1.0)

        ! Primer numero aleatorio
        I1 = mod(a0 * s + c0, m0)

        ! Evitar x = 0
        if (I1 == 0) I1 = 1

        ! Segundo numero aleatorio
        I2 = mod(a0 * I1 + c0, m0)

        ! Evitar y = 0
        if (I2 == 0) I2 = 1

        ! Guardamos la siguiente semilla
        In2 = I2

        ! Normalizacion
        x = real(I1) / real(m0)
        y = real(I2) / real(m0)

        ! Box-Muller
        no1 = - sqrt(-2.0 * log(x)) * cos(2.0 * pi * y)

        no2 = - sqrt(-2.0 * log(x)) * sin(2.0 * pi * y)

    end subroutine generar_normal


    subroutine generar_normal_txt_congr(pasos, semilla, a, c, m, nombre_archivo)

        implicit none

        integer(kind=8), intent(in) :: a, c, m
        integer(kind=8), intent(inout) :: semilla
        integer, intent(in) :: pasos
        character(len=*), intent(in) :: nombre_archivo

        integer(kind=8) :: aux
        integer :: i
        real :: n1, n2

        open(unit=2, file=nombre_archivo, status="replace")
        do i = 1, pasos

            call generar_normal(semilla, a, c, m, n1, n2, aux)

            ! La siguiente semilla
            semilla = aux

            write(2,*) n1, n2
        end do
        close(2)

    end subroutine generar_normal_txt_congr

    subroutine estudio_rand_f90(iteraciones, nombre_archivo)
        character(len=*), intent(in) :: nombre_archivo
        integer, intent(in) :: iteraciones
        real :: x,y
        integer :: i

        call random_number(x)

        open(unit=10, file=nombre_archivo, status="replace")

         do i=1, iteraciones    
            call random_number(y)
            write(10,*) x, y

            x=y
        
         end do
    
    end subroutine estudio_rand_f90

    subroutine generar_normal_txt_f90(pasos, nombre_archivo)

        implicit none

        integer, intent(in) :: pasos
        character(len=*), intent(in) :: nombre_archivo
        integer :: i
        real :: q1, q2

        open(unit=13, file=nombre_archivo, status="replace")
        do i = 1, pasos

            call generar_normal_f90(q1,q2)

            write(13,*) q1, q2
        end do
        close(13)
    
    end subroutine generar_normal_txt_f90

    subroutine generar_normal_f90(q1,q2)

        implicit none
        
        real, intent(out) :: q1,q2
        real :: pi, x, y

        pi = 4.0 * atan(1.0)

        ! Normalizacion
        call random_number (x)
        call random_number (y)

        if (x==0) x=1
        if (y==0) y=1

        ! Box-Muller
        q1 = - sqrt(-2.0 * log(x)) * cos(2.0 * pi * y)

        q2 = - sqrt(-2.0 * log(x)) * sin(2.0 * pi * y)

    end subroutine generar_normal_f90

    subroutine estudio_rand_congr(s, a, c, m, iteraciones, nombre_archivo)

        implicit none

        integer(kind=8), intent(in) :: s, a, c, m
        integer, intent(in) :: iteraciones
        character(len=*), intent(in) :: nombre_archivo
        integer (kind=8) :: I1
        real :: x, y
        integer :: i


        open(unit=11, file=nombre_archivo, status="replace")


        I1 = mod(a * s + c, m)
        x=real(I1)/real(m)

        do i=1, iteraciones 
            I1=mod(a * I1 + c,m)
            y=real(I1)/real(m)

            write(11,*) x,y

            x=y
            
        end do
        
    end subroutine estudio_rand_congr

end module generador_normal

!========================== euler_maruyama(pasos,h,xo,po,m,k,beta_inv,nu) ===========================================
!========================== Runge_Kutta_2 (pasos,h,xo,po,m,k,beta_inv,nu) ===========================================
!========================== verlet_exp_est(pasos,h,xo,po,m,k,beta_inv,nu) ===========================================
module algoritmos_estocastios
    use generador_normal
    contains

    subroutine euler_maruyama_01(pasos,h,xo,po,m,k,beta_inv,nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        real :: x, p, n1, n2, aux, Ki, V, E
        integer :: i

        x = xo
        p = po

        open(unit=3, file="euler_maruyama.txt", status="replace")
        Ki=0.5*p*p/m
        V=0.5*k*x*x
        E=Ki+V
        write(3,*) 0.0, x, p, Ki, V, E
        do i = 1, pasos/2

            call generar_normal_f90(n1, n2)
            
            ! Actualizacion de las variables usando Euler-Maruyama
            aux = x
            x = x + (p/m)*h 
            p = p - k*aux*h - nu*p*h + sqrt(2.0*nu*m*beta_inv*h)*n1

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(3,*) (2*i-1)*h,x, p, Ki, V, E

            aux = x
            x = x + (p/m)*h 
            p = p - k*aux*h - nu*p*h + sqrt(2.0*nu*m*beta_inv*h)*n2

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(3,*) (2*i)*h, x, p, Ki, V, E

        end do
        close(3)

    end subroutine euler_maruyama_01

    subroutine Runge_Kutta_2_01 (pasos, h, xo, po, m, k, beta_inv, nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        real :: x, p, n1, n2, fx1, fx2, fg1, fg2, z_det, Ki, V, E
        integer :: i

        x = xo
        p = po

        open(unit=1001, file="runge_kutta_2.txt", status="replace")
        Ki=0.5*p*p/m
        V=0.5*k*x*x
        E=Ki+V
        write(1001,*) 0.0,x, p, Ki, V, E
        z_det=sqrt(2.0*nu*m*beta_inv*h)
        do i = 1, pasos/2

            call generar_normal_f90(n1, n2)

            fx1 = (p+z_det*n1)/m
            fg1 = -k*x - nu*(p+z_det*n1)

            fx2 = (p + h*fg1)/m
            fg2 = -k*(x + h*fx1) - nu*(p + h*fg1)

            x = x + 0.5*h*(fx1 + fx2)
            p = p + 0.5*h*(fg1 + fg2) + z_det*n1

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(1001,*) (2*i-1)*h, x, p, Ki, V, E

            fx1 = (p+z_det*n2)/m
            fg1 = -k*x - nu*(p+z_det*n2)

            fx2 = (p + h*fg1)/m
            fg2 = -k*(x + h*fx1) - nu*(p + h*fg1)

            x = x + 0.5*h*(fx1 + fx2)
            p = p + 0.5*h*(fg1 + fg2) + z_det*n2

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(1001,*) (2*i)*h, x, p, Ki, V, E

        end do

        close(1001)

    end subroutine Runge_Kutta_2_01

    subroutine verlet_exp_est_01(pasos, h, xo, po, m, k, beta_inv, nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        real :: x, p, n1, n2, z_det, verlet_a, verlet_b, aux, Ki, V, E
        integer :: i

        x = xo
        p = po

        open(unit=5, file="verlet_exp_est.txt", status="replace")
        Ki=0.5*p*p/m
        V=0.5*k*x*x
        E=Ki+V
        write(5,*) 0.0,x, p, Ki, V, E
        z_det=sqrt(2.0*nu*m*beta_inv*h)
        verlet_a = (1-(nu*h)/2)/(1+(nu*h)/2)
        verlet_b = 1/(1+(nu*h)/2)

        do i = 1, pasos/2

            call generar_normal_f90(n1, n2)

            aux = x
            x = x + verlet_b*h*p/m - verlet_b*h*h*k*x/(2.0*m) + verlet_b*h*z_det*n1/(2.0*m)
            p = verlet_a*p - 0.5*h*(verlet_a*k*aux + k*x) + verlet_b*z_det*n1

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(5,*) (2*i-1)*h, x, p, Ki, V, E

            aux = x
            x = x + verlet_b*h*p/m - verlet_b*h*h*k*x/(2.0*m) + verlet_b*h*z_det*n2/(2.0*m)
            p = verlet_a*p - 0.5*h*(verlet_a*k*aux + k*x) + verlet_b*z_det*n2

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(5,*) (2*i)*h, x, p, Ki, V, E

        end do
        close(5)

    end subroutine verlet_exp_est_01

    subroutine termalizacion_EM(pasos, h, xo, po, m, k, beta_inv, nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        real :: x, p, n1, n2, aux, Ki, V, Ki_m, V_m
        integer :: i

        x = xo
        p = po

        open(unit=15, file="termalizacion_EM.txt", status="replace")
        Ki=0.5*p*p/m
        V=0.5*k*x*x
        Ki_m = Ki
        V_m = V
        write(15,*) 0.0, x, p, Ki_m, V_m
        do i = 1, pasos/2
            call generar_normal_f90(n1,n2)
            aux = x
            x = x + (p/m)*h 
            p = p - k*aux*h - nu*p*h + sqrt(2.0*nu*m*beta_inv*h)*n1

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            Ki_m = Ki_m + Ki
            V_m=V_m + V

            write(15,*) (2*i-1)*h,x, p, Ki_m/(2*i), V_m/(2*i)

            aux = x
            x = x + (p/m)*h 
            p = p - k*aux*h - nu*p*h + sqrt(2.0*nu*m*beta_inv*h)*n2

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            Ki_m = Ki_m + Ki
            V_m = V_m + V

            write(15,*) (2*i)*h,x, p, Ki_m/((2*i)+1), V_m/((2*i)+1)
        end do
        close (15)


        
    end subroutine termalizacion_EM

    subroutine termalizacion_RK(pasos, h, xo, po, m, k, beta_inv, nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        real :: x, p, n1, n2, fx1, fx2, fg1, fg2, z_det, Ki, V, Ki_m, V_m
        integer :: i

        x = xo
        p = po

        open(unit=16, file="term_RK.txt", status="replace")
        Ki=0.5*p*p/m
        V=0.5*k*x*x
        Ki_m=Ki
        V_m=V
        write(16,*) 0.0,x, p, Ki_m, V_m
        z_det=sqrt(2.0*nu*m*beta_inv*h)
        do i = 1, pasos/2

            call generar_normal_f90(n1, n2)

            fx1 = (p+z_det*n1)/m
            fg1 = -k*x - nu*(p+z_det*n1)

            fx2 = (p + h*fg1)/m
            fg2 = -k*(x + h*fx1) - nu*(p + h*fg1)

            x = x + 0.5*h*(fx1 + fx2)
            p = p + 0.5*h*(fg1 + fg2) + z_det*n1

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            Ki_m = Ki_m + Ki
            V_m = V_m + V
            write(16,*) (2*i-1)*h, x, p, Ki_m/(2*i), V_m/(2*i)

            fx1 = (p+z_det*n2)/m
            fg1 = -k*x - nu*(p+z_det*n2)

            fx2 = (p + h*fg1)/m
            fg2 = -k*(x + h*fx1) - nu*(p + h*fg1)

            x = x + 0.5*h*(fx1 + fx2)
            p = p + 0.5*h*(fg1 + fg2) + z_det*n2

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            Ki_m=Ki_m + Ki
            V_m = V_m + V
            write(16,*) (2*i)*h, x, p, Ki_m/((2*i)+1), V_m/((2*i)+1)

        end do

        close(16)
    
    end subroutine termalizacion_RK

    subroutine termalizacion_VE(pasos, h, xo, po, m, k, beta_inv, nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        real :: x, p, n1, n2, z_det, verlet_a, verlet_b, aux, Ki, V, Ki_m, V_m
        integer :: i

        x = xo
        p = po

        open(unit=17, file="term_VE.txt", status="replace")
        Ki=0.5*p*p/m
        V=0.5*k*x*x
        Ki_m=Ki
        V_m=V
        write(17,*) 0.0,x, p, Ki_m, V_m
        z_det=sqrt(2.0*nu*m*beta_inv*h)
        verlet_a = (1-(nu*h)/2)/(1+(nu*h)/2)
        verlet_b = 1/(1+(nu*h)/2)

        do i = 1, pasos/2

            call generar_normal_f90(n1, n2)

            aux = x
            x = x + verlet_b*h*p/m - verlet_b*h*h*k*x/(2.0*m) + verlet_b*h*z_det*n1/(2.0*m)
            p = verlet_a*p - 0.5*h*(verlet_a*k*aux + k*x) + verlet_b*z_det*n1

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            Ki_m=Ki_m + Ki
            V_m = V_m + V
            write(17,*) (2*i-1)*h, x, p, Ki_m/(2*i), V_m/(2*i)

            aux = x
            x = x + verlet_b*h*p/m - verlet_b*h*h*k*x/(2.0*m) + verlet_b*h*z_det*n2/(2.0*m)
            p = verlet_a*p - 0.5*h*(verlet_a*k*aux + k*x) + verlet_b*z_det*n2

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            Ki_m=Ki_m + Ki
            V_m = V_m + V
            write(17,*) (2*i)*h, x, p, Ki_m/((2*i)+1), V_m/((2*i)+1)

        end do
        close(17)

    end subroutine termalizacion_VE

end module algoritmos_estocastios

module objetivo_2
    use generador_normal
    contains

    subroutine verlet_doblepozo(pasos,h,xo,po,m,A,beta_inv,eta,F)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, A, beta_inv, eta, F
        real(kind=8) :: Ki_m, V_m, Ki, V
        real :: x, p, n1, n2, z_det, verlet_a, verlet_b, aux, t_min
        integer :: i, Contador, Contador2
        
        Contador2=0
        Contador=0
        x = xo
        p = po
        Ki=0.5*p*p/m
        V=A*(x**2-1)**2
        Ki_m=Ki
        V_m=V
        if (x>0) then
            Contador=Contador+1
        end if

        if (eta < 0.1 .and. abs(A-0.5) < 1.0e-6) then
            t_min = 150.0
        else if (eta < 0.1 .and. abs(A-1.0) < 1.0e-6) then
            t_min = 300.0        ! eta pequeño: muchos recruces, filtro grande
        else if (eta < 1.0 .and. abs(A-0.5) < 1.0e-6) then
            t_min = 20.0         ! eta intermedio
        else if (eta < 1.0 .and. abs(A-1.0) < 1.0e-6) then
            t_min = 50.0
        else
            t_min = 50.0
        end if

        open(unit=176, file="verlet_doblepozo.txt", status="replace")
        open(unit=177, file="t_estancia.txt", status="replace")

        write(176,*) 0.0,x,p, Ki_m, V_m, Contador
        z_det=sqrt(2.0*eta*m*beta_inv*h)
        verlet_a = (1-(eta*h)/2)/(1+(eta*h)/2)
        verlet_b = 1/(1+(eta*h)/2)

        do i = 1, pasos/2

            call generar_normal_f90(n1,n2)

            aux = x
            x = x + verlet_b*h*p/m + verlet_b*h*h*(-4*A*x*(x**2-1)+F)/(2.0*m) + verlet_b*h*z_det*n1/(2.0*m)
            p = verlet_a*p + h*(verlet_a*(-4*A*aux*(aux**2-1)+F) + (-4*A*x*(x**2-1)+F))/2 + verlet_b*z_det*n1

            Ki=0.5*p*p/m
            V=A*(x**2-1)**2
            Ki_m=Ki_m + Ki
            V_m = V_m + V
            if (x>0.0) then
                Contador=Contador+1
            end if
            write(176,*) (2*i-1)*h,x,p, Ki_m/(2*i), V_m/(2*i), real(Contador)/(real(2*i))

            if(x*aux>0.0) then
                Contador2=Contador2+1
            else
                if (contador2*h>t_min) write(177,*) Contador2*h
                Contador2 = 0
            end if

            aux = x
            x = x + verlet_b*h*p/m + verlet_b*h*h*(-4*A*x*(x**2-1)+F)/(2*m) + verlet_b*h*z_det*n2/(2.0*m)
            p = verlet_a*p + h*(verlet_a*(-4*A*aux*(aux**2-1)+F) + (-4*A*x*(x**2-1)+F))/2 + verlet_b*z_det*n2

            Ki=0.5*p*p/m
            V=A*(x**2-1)**2
            Ki_m=Ki_m + Ki
            V_m = V_m + V
            if (x>0.0) then
                Contador=Contador+1
            end if
            write(176,*) (2*i)*h,x,p, Ki_m/((2*i)+1), V_m/((2*i)+1), real(Contador)/(real((2*i)+1))

            if(x*aux>0.0) then
                Contador2=Contador2+1
            else
                if (contador2*h>t_min) write(177,*) Contador2*h
                Contador2 = 0
            end if

        end do

        close(176)
        close(177)
    end subroutine verlet_doblepozo

    subroutine verlet_doblepozo_testancia(pasos,h,xo,po,m,A,beta_inv,eta,F)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, A, beta_inv, eta, F
        real :: x, p, n1, n2, z_det, verlet_a, verlet_b, aux, t_min
        integer :: i, Contador
        
        Contador=0
        x = xo
        p = po

        if (eta < 0.1 .and. abs(A-0.5) < 1.0e-6) then
            t_min = 150.0
        else if (eta < 0.1 .and. abs(A-1.0) < 1.0e-6) then
            t_min = 300.0        ! eta pequeño: muchos recruces, filtro grande
        else if (eta < 1.0 .and. abs(A-0.5) < 1.0e-6) then
            t_min = 20.0         ! eta intermedio
        else if (eta < 1.0 .and. abs(A-1.0) < 1.0e-6) then
            t_min = 50.0
        else
            t_min = 50.0
        end if

        open(unit=200, file="t_estancia_only.txt", status="replace")
        z_det=sqrt(2.0*eta*m*beta_inv*h)
        verlet_a = (1-(eta*h)/2)/(1+(eta*h)/2)
        verlet_b = 1/(1+(eta*h)/2)

        do i = 1, pasos/2

            call generar_normal_f90(n1,n2)

            aux = x
            x = x + verlet_b*h*p/m + verlet_b*h*h*(-4*A*x*(x**2-1)+F)/(2.0*m) + verlet_b*h*z_det*n1/(2.0*m)
            p = verlet_a*p + h*(verlet_a*(-4*A*aux*(aux**2-1)+F) + (-4*A*x*(x**2-1)+F))/2 + verlet_b*z_det*n1

            if(x*aux>0.0) then
                Contador=Contador+1
            else
                if (contador*h>t_min) write(200,*) Contador*h
                Contador = 0
            end if

            aux = x
            x = x + verlet_b*h*p/m + verlet_b*h*h*(-4*A*x*(x**2-1)+F)/(2*m) + verlet_b*h*z_det*n2/(2.0*m)
            p = verlet_a*p + h*(verlet_a*(-4*A*aux*(aux**2-1)+F) + (-4*A*x*(x**2-1)+F))/2 + verlet_b*z_det*n2

            if(x*aux>0.0) then
                Contador=Contador+1
            else
                if (contador*h>t_min) write(200,*) Contador*h
                Contador = 0
            end if

        end do

        close(200)
    end subroutine verlet_doblepozo_testancia

    subroutine verlet_dp_fuerza(pasos,h,xo,po,m,A,beta_inv,eta,F)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, A, beta_inv, eta, F
        real :: x, p, n1, n2, z_det, verlet_a, verlet_b, aux
        integer :: i, Contador
        
        Contador=0
        x = xo
        p = po
        if (x>0) then
            Contador=Contador+1
        end if

        z_det=sqrt(2.0*eta*m*beta_inv*h)
        verlet_a = (1-(eta*h)/2)/(1+(eta*h)/2)
        verlet_b = 1/(1+(eta*h)/2)

        do i = 1, pasos/2

            call generar_normal_f90(n1,n2)

            aux = x
            x = x + verlet_b*h*p/m + verlet_b*h*h*(-4*A*x*(x**2-1)+F)/(2.0*m) + verlet_b*h*z_det*n1/(2.0*m)
            p = verlet_a*p + h*(verlet_a*(-4*A*aux*(aux**2-1)+F) + (-4*A*x*(x**2-1)+F))/2 + verlet_b*z_det*n1

            if (x>0.0) then
                Contador=Contador+1
            end if

            aux = x
            x = x + verlet_b*h*p/m + verlet_b*h*h*(-4*A*x*(x**2-1)+F)/(2*m) + verlet_b*h*z_det*n2/(2.0*m)
            p = verlet_a*p + h*(verlet_a*(-4*A*aux*(aux**2-1)+F) + (-4*A*x*(x**2-1)+F))/2 + verlet_b*z_det*n2

            if (x>0.0) then
                Contador=Contador+1
            end if

        end do

        print *, "Fraccion de tiempo en el pozo 2: ", real(Contador)/real(pasos)

    end subroutine verlet_dp_fuerza

    

end module objetivo_2

module objetivo_3
    use generador_normal
    implicit none
    integer, parameter :: dp = kind(1.0d0)   ! precision doble (main la hereda con "use objetivo_3")
    contains

    subroutine fuerzas(r,f,ue,l0,k)
        real(dp), intent(in) :: r(:, :), l0, k
        real(dp), intent(out) :: f(:, :), ue
        real(dp) :: d(3), l
        integer :: i,n
        n=size(r,2)
        f=0.0_dp
        ue=0.0_dp
        do i=1,n-1
            d = r(:,i+1) - r(:,i)
            l=sqrt(sum(d**2))
            ue=ue+0.5_dp*k*(l-l0)**2
            f(:,i+1)=f(:,i+1)-k*(l-l0)*(d/l)
            f(:,i)=f(:,i)+k*(l-l0)*(d/l)
        end do
    end subroutine fuerzas

    subroutine polimero(pasos, h, r, p, m, k, beta_inv, eta, l0)
        integer, intent(in) :: pasos
        real(dp), intent(in) :: h, m, k, beta_inv, eta, l0
        real(dp), intent(inout) :: r(:, :), p(:, :)
        real(dp) :: f(size(r,1), size(r,2)), aux(size(r,1),size(r,2))
        real(dp) :: n1(size(r,1),size(r,2)), n2(size(r,1),size(r,2))
        real(dp) :: z_det, verlet_a, verlet_b
        real :: gaus1, gaus2
        real(dp) :: ue, Ki_m, V_m
        integer :: i,j,ka

        open(unit=180, file="Termalizacion_polimero.txt", status="replace")

        z_det=sqrt(2.0_dp*eta*m*beta_inv*h)
        verlet_a = (1.0_dp-(eta*h)/2.0_dp)/(1.0_dp+(eta*h)/2.0_dp)
        verlet_b = 1.0_dp/(1.0_dp+(eta*h)/2.0_dp)
        call fuerzas(r,f,ue,l0,k)

        ki_m=0.5_dp*sum(p**2)/m
        V_m=ue

        write(180,*) 0.0_dp,Ki_m, V_m


        do i = 1, pasos/2

            do ka=1,size(r,2)
                do j=1, size(r,1)
                    call generar_normal_f90(gaus1,gaus2)
                    n1(j,ka)=gaus1
                    n2(j,ka)=gaus2
                end do 
            end do 

            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n1/(2.0_dp*m)
            call fuerzas(r,f,ue,l0,k)
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n1
            
            Ki_m=Ki_m+0.5_dp*sum(p**2)/m
            V_m=V_m+ue

            if (mod(i,50) == 0) write(180,*) (2*i-1)*h, Ki_m/(2*i), V_m/(2*i)

            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n2/(2.0_dp*m)
            call fuerzas(r,f,ue,l0,k)
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n2

            Ki_m=Ki_m+0.5_dp*sum(p**2)/m
            V_m=V_m+ue

            if (mod(i,50) == 0) write(180,*) 2*i*h, Ki_m/(2*i+1), V_m/(2*i+1)

        end do

        close(180)

    end subroutine polimero
    
    subroutine radio_giro (pasos, h, r, p, m, k, beta_inv, eta, l0)
        integer, intent(in) :: pasos
        real(dp), intent(in) :: h, m, k, beta_inv, eta, l0
        real(dp), intent(inout) :: r(:, :), p(:, :)
        real(dp) :: f(size(r,1), size(r,2)), aux(size(r,1),size(r,2))
        real(dp) :: n1(size(r,1),size(r,2)), n2(size(r,1),size(r,2))
        real(dp) :: z_det, verlet_a, verlet_b
        real :: gaus1, gaus2
        real(dp) :: rcm(3), rg2_m, rg_m, ue, rg2, rg
        integer :: i,j,ka

        !rg2_m, rg_m, raiz(rg2_m), rg_m*2
        open(unit=181, file="Radio_giro.txt", status="replace")

        z_det=sqrt(2.0_dp*eta*m*beta_inv*h)
        verlet_a = (1.0_dp-(eta*h)/2.0_dp)/(1.0_dp+(eta*h)/2.0_dp)
        verlet_b = 1.0_dp/(1.0_dp+(eta*h)/2.0_dp)
        call fuerzas(r,f,ue,l0,k)

        rcm = sum(r, dim=2) / size(r,2)
        !Calculo radio de giro al cuadrado del polimero en la configuracion inicial
        rg2 = 0.0_dp
        rg= 0.0_dp
        do i = 1, size(r,2)
            rg2 = rg2 + sum((r(:,i) - rcm)**2)
        end do
        rg2 = rg2 / size(r,2)

        rg=sqrt(rg2)

        rg2_m=rg2
        rg_m=rg

        write(181,*) 0.0_dp, rg2_m, rg_m, sqrt(rg2_m), (rg_m)**2


        do i = 1, pasos/2

            do ka=1,size(r,2)
                do j=1, size(r,1)
                    call generar_normal_f90(gaus1,gaus2)
                    n1(j,ka)=gaus1
                    n2(j,ka)=gaus2
                end do 
            end do 

            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n1/(2.0_dp*m)
            call fuerzas(r,f,ue,l0,k)
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n1
            
            rcm = sum(r, dim=2) / size(r,2)
            rg2 = 0.0_dp
            do j = 1, size(r,2)
                rg2 = rg2 + sum((r(:,j) - rcm)**2)
            end do
            rg2 = rg2 / size(r,2)

            rg = sqrt(rg2)

            rg2_m=rg2_m + rg2
            rg_m=rg_m + rg
            
            !Guardamos en el archivo de salida el tiempo y el radio de giro promedio cada 50 pasos
            if (mod(i,50) == 0) write(181,*) (2*i-1)*h, rg2_m/(2*i), rg_m/(2*i), sqrt(rg2_m/(2*i)), (rg_m/(2*i))**2

            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n2/(2.0_dp*m)
            call fuerzas(r,f,ue,l0,k)
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n2

            rcm = sum(r, dim=2) / size(r,2)
            rg2 = 0.0_dp
            do j = 1, size(r,2)
                rg2 = rg2 + sum((r(:,j) - rcm)**2)
            end do
            rg2 = rg2 / size(r,2)

            rg = sqrt(rg2)

            rg2_m=rg2_m + rg2

            rg_m = rg_m + rg

            if (mod(i,50) == 0) write(181,*) 2*i*h, rg2_m/(2*i+1), rg_m/(2*i+1), sqrt(rg2_m/(2*i+1)), (rg_m/(2*i+1))**2

        end do

        close(181)

    end subroutine radio_giro

    subroutine dist_extremos (pasos, h, r, p, m, k, beta_inv, eta, l0)
        integer, intent(in) :: pasos
        real(dp), intent(in) :: h, m, k, beta_inv, eta, l0
        real(dp), intent(inout) :: r(:, :), p(:, :)
        real(dp) :: f(size(r,1), size(r,2)), aux(size(r,1),size(r,2))
        real(dp) :: n1(size(r,1),size(r,2)), n2(size(r,1),size(r,2))
        real(dp) :: z_det, verlet_a, verlet_b
        real :: gaus1, gaus2
        real(dp) :: ue, Ree, Ree_m
        integer :: i,j,ka

        !rg2_m, rg_m, raiz(rg2_m), rg_m*2
        open(unit=197, file="Radio_giro.txt", status="replace")

        z_det=sqrt(2.0_dp*eta*m*beta_inv*h)
        verlet_a = (1.0_dp-(eta*h)/2.0_dp)/(1.0_dp+(eta*h)/2.0_dp)
        verlet_b = 1.0_dp/(1.0_dp+(eta*h)/2.0_dp)
        call fuerzas(r,f,ue,l0,k)

        Ree= sqrt(sum((r(:,1)-r(:,size(r,2)))**2))
        Ree_m=Ree

        write(197,*) 0.0_dp, Ree


        do i = 1, pasos/2

            do ka=1,size(r,2)
                do j=1, size(r,1)
                    call generar_normal_f90(gaus1,gaus2)
                    n1(j,ka)=gaus1
                    n2(j,ka)=gaus2
                end do 
            end do 

            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n1/(2.0_dp*m)
            call fuerzas(r,f,ue,l0,k)
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n1
            
            Ree= sqrt(sum((r(:,1)-r(:,size(r,2)))**2))
            Ree_m=Ree_m + Ree
            if (mod(i,50) == 0) write(197,*) (2*i-1)*h, Ree_m/(2*i)

            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n2/(2.0_dp*m)
            call fuerzas(r,f,ue,l0,k)
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n2

            Ree= sqrt(sum((r(:,1)-r(:,size(r,2)))**2))
            Ree_m=Ree_m + Ree
            if (mod(i,50) == 0) write(197,*) (2*i)*h, Ree_m/((2*i)+1)

        end do

        close(197)

    end subroutine dist_extremos

    subroutine dist_extremos_F (pasos, h, r, p, m, k, beta_inv, eta, l0, Fuerza)
        integer, intent(in) :: pasos
        real(dp), intent(in) :: h, m, k, beta_inv, eta, l0
        real(dp), intent(inout) :: r(:, :), p(:, :)
        real(dp), intent(in) :: Fuerza(size(r,1)) 
        real(dp) :: r0(size(r,1))
        real(dp) :: f(size(r,1), size(r,2)), aux(size(r,1),size(r,2))
        real(dp) :: n1(size(r,1),size(r,2)), n2(size(r,1),size(r,2))
        real(dp) :: z_det, verlet_a, verlet_b
        real :: gaus1, gaus2
        real(dp) :: ue, Ree, Ree_m
        integer :: i,j,ka

        !rg2_m, rg_m, raiz(rg2_m), rg_m*2
        open(unit=198, file="Radio_giro_F.txt", status="replace")
        
        r0=r(:,1)
        z_det=sqrt(2.0_dp*eta*m*beta_inv*h)
        verlet_a = (1.0_dp-(eta*h)/2.0_dp)/(1.0_dp+(eta*h)/2.0_dp)
        verlet_b = 1.0_dp/(1.0_dp+(eta*h)/2.0_dp)
        call fuerzas(r,f,ue,l0,k)
        f(:,size(r,2))=f(:,size(r,2)) + Fuerza

        Ree = sum((r(:,size(r,2))-r(:,1))*Fuerza) / sqrt(sum(Fuerza**2)) 
        Ree_m=Ree
        write(198,*) 0.0_dp, Ree

        do i = 1, pasos/2
            do ka=1,size(r,2)
                do j=1, size(r,1)
                    call generar_normal_f90(gaus1,gaus2)
                    n1(j,ka)=gaus1
                    n2(j,ka)=gaus2
                end do 
            end do 

            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n1/(2.0_dp*m)
            r(:,1)=r0
            call fuerzas(r,f,ue,l0,k)
            f(:,size(r,2)) = f(:,size(r,2)) + Fuerza
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n1
            
            Ree = sum((r(:,size(r,2))-r(:,1))*Fuerza) / sqrt(sum(Fuerza**2)) 
            Ree_m=Ree_m + Ree
            if (mod(i,50) == 0) write(198,*) (2*i-1)*h, Ree_m/(2*i)

            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n2/(2.0_dp*m)
            r(:,1)=r0
            call fuerzas(r,f,ue,l0,k)
            f(:,size(r,2))=f(:,size(r,2)) + Fuerza
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n2

            Ree = sum((r(:,size(r,2))-r(:,1))*Fuerza) / sqrt(sum(Fuerza**2)) 
            Ree_m=Ree_m + Ree
            if (mod(i,50) == 0) write(198,*) (2*i)*h, Ree_m/((2*i)+1)

        end do

        close(198)

    end subroutine dist_extremos_F

    subroutine VMD(pasos, h, r, p, m, k, beta_inv, eta, l0, dataframe)
        integer, intent(in) :: pasos, dataframe
        real(dp), intent(in) :: h, m, k, beta_inv, eta, l0
        real(dp), intent(inout) :: r(:, :), p(:, :)
        real(dp) :: f(size(r,1), size(r,2)), aux(size(r,1),size(r,2))
        real(dp) :: n1(size(r,1),size(r,2)), n2(size(r,1),size(r,2))
        real(dp) :: z_det, verlet_a, verlet_b, t, ue
        real :: gaus1, gaus2
        integer :: i,j,ka

        !Aquí guardaremos la energía cinética y potencial promedio del polimero a lo largo de la simulación
        open(unit=190, file="Frames_VMD.xyz", status="replace")

        t=0.0_dp

        z_det=sqrt(2.0_dp*eta*m*beta_inv*h)
        verlet_a = (1.0_dp-(eta*h)/2.0_dp)/(1.0_dp+(eta*h)/2.0_dp)
        verlet_b = 1.0_dp/(1.0_dp+(eta*h)/2.0_dp)
        call fuerzas(r,f,ue,l0,k)
        !Energía cinética y potencial inicial del polimero

        write(190,'(I0)') size(r,2)
        write(190,'(A,F12.3,A,F8.2,A,F6.2,A,F6.2)') 't=', t, '  ke=', k, '  l0=', l0, '  eta=',eta
        do j=1, size(r,2)
            write(190,'(A,3F12.5)') 'C ', r(1,j), r(2,j), r(3,j)
        end do

        !Bucle principal de la simulación, donde se actualizan las posiciones y velocidades del polimero usando el algoritmo de Verlet estocástico
        !Hace pasos/2 iteraciones, y en cada iteración se generan dos pasos de tiempo, uno con n1 y otro con n2
        do i = 1, pasos/2

            !Generamos dos matrices de números aleatorios gaussianos n1 y n2, que se usarán para actualizar las posiciones y velocidades del polimero
            do ka=1,size(r,2)
                do j=1, size(r,1)
                    call generar_normal_f90(gaus1,gaus2)
                    n1(j,ka)=gaus1
                    n2(j,ka)=gaus2
                end do 
            end do 
            
            !Actualizamos las posiciones y velocidades del polimero usando el algoritmo de Verlet estocástico con n1
            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n1/(2.0_dp*m)
            call fuerzas(r,f,ue,l0,k)
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n1
            
            !Calculamos la energía cinética y potencial del polimero después de este paso, y las acumulamos en Ki_m y V_m para obtener el promedio
            t = t + h
            if (mod((2*i-1),dataframe) == 0)  then
                write(190,'(I0)') size(r,2)
                write(190,'(A,F12.3,A,F8.2,A,F6.2,A,F6.2)') 't=', t, '  ke=', k, '  l0=', l0, '  eta=',eta
                do j=1, size(r,2)
                    write(190,'(A,3F12.5)') 'C ', r(1,j), r(2,j), r(3,j)
                end do
            end if

            !Repetimos el mismo proceso con n2, generando otro paso de tiempo y actualizando las posiciones y velocidades del polimero
            aux = f
            r = r + verlet_b*h*p/m + verlet_b*h*h*f/(2.0_dp*m) + verlet_b*h*z_det*n2/(2.0_dp*m)
            call fuerzas(r,f,ue,l0,k)
            p = verlet_a*p + h*(verlet_a*aux + f)/2.0_dp + verlet_b*z_det*n2

            t = t + h
            if (mod((2*i),dataframe) == 0)  then
                write(190,'(I0)') size(r,2)
                write(190,'(A,F12.3,A,F8.2,A,F6.2,A,F6.2)') 't=', t, '  ke=', k, '  l0=', l0, '  eta=',eta
                do j=1, size(r,2)
                    write(190,'(A,3F12.5)') 'C ', r(1,j), r(2,j), r(3,j)
                end do
                
            end if

        end do

        close(190)

    end subroutine VMD

end module objetivo_3


program main
    use objetivo_3, only: dp, fuerzas
    implicit none

    integer, parameter :: tamanos(*) = [4, 6, 8, 10, 12, 16, 20, &
                                         24, 28, 32, 40, 48, 64]
    real(dp), parameter :: h = 0.001_dp, m = 1.0_dp
    real(dp), parameter :: ke = 1000.0_dp, kBT = 1.0_dp
    real(dp), parameter :: eta = 1.0_dp, l0 = 1.0_dp
    real(dp) :: t_termal, t_medicion
    integer, parameter :: cada_muestra = 500
    integer, parameter :: cada_guardado = 500000
    real(dp), parameter :: pi = acos(-1.0_dp)

    integer :: j, pasos_termal, pasos_medicion
    real(dp) :: a, b, amplitud_ruido
    character(len=240) :: archivo
    logical :: existe

    if (mod(cada_guardado, cada_muestra) /= 0) then
        error stop 'Frecuencias de muestreo incompatibles'
    end if

    a = (1.0_dp - eta*h/2.0_dp) / (1.0_dp + eta*h/2.0_dp)
    b = 1.0_dp / (1.0_dp + eta*h/2.0_dp)
    amplitud_ruido = sqrt(2.0_dp*eta*m*kBT*h)

    ! Verificar todas las salidas antes de comenzar un calculo largo.
    do j = 1, size(tamanos)
        call ruta_salida(tamanos(j), archivo)
        inquire(file=trim(archivo), exist=existe)
        if (existe) then
            print *, 'Archivo ya existente: ', trim(archivo)
            error stop 'No se sobreescribiran los resultados'
        end if
    end do

    call random_seed()
    print *, 'RADIO DE GIRO Y DISTANCIA ENTRE EXTREMOS - ke=', ke
    print *, 'N:', tamanos

    do j = 1, size(tamanos)
        ! Dedicar mas tiempo a las cadenas largas.
        if (tamanos(j) == 64) then
            t_termal = 40000.0_dp
            t_medicion = 200000.0_dp
        else if (tamanos(j) >= 40) then
            t_termal = 20000.0_dp
            t_medicion = 150000.0_dp
        else
            t_termal = 10000.0_dp
            t_medicion = 100000.0_dp
        end if

        pasos_termal = nint(t_termal/h)
        pasos_medicion = nint(t_medicion/h)
        if (mod(pasos_medicion, cada_guardado) /= 0) then
            error stop 'Tiempo de medicion incompatible con guardado'
        end if
        print *, 'N =', tamanos(j), ' termal =', t_termal, &
                 ' medicion =', t_medicion
        call simular_cadena(tamanos(j))
    end do
    print *, 'Barrido terminado'

contains

    subroutine ruta_salida(n, nombre)
        integer, intent(in) :: n
        character(len=*), intent(out) :: nombre
        character(len=16) :: texto_n
        write(texto_n,'(I0)') n
        nombre = 'Py/plots/OBJETIVO_3/Datos/k1000/RG/' // &
                 'Rg_Ree_N' // trim(texto_n) // '_ke1000.txt'
    end subroutine ruta_salida

    subroutine normales(z1,z2)
        real(dp), intent(out) :: z1,z2
        real(dp) :: u1,u2,rho,theta
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1,tiny(1.0_dp))
        rho = sqrt(-2.0_dp*log(u1))
        theta = 2.0_dp*pi*u2
        z1 = rho*cos(theta)
        z2 = rho*sin(theta)
    end subroutine normales

    subroutine inicializar(n,r,p)
        integer, intent(in) :: n
        real(dp), intent(out) :: r(3,n),p(3,n)
        real(dp) :: z1,z2,z3,z4,orientacion(3),longitud
        integer :: i,j
        r(:,1) = 0.0_dp
        do i = 2,n
            call normales(z1,z2)
            call normales(z3,z4)
            orientacion = [z1,z2,z3]
            orientacion = orientacion/sqrt(sum(orientacion**2))
            longitud = l0 + sqrt(kBT/ke)*z4
            do while (longitud <= 0.0_dp)
                call normales(z3,z4)
                longitud = l0 + sqrt(kBT/ke)*z4
            end do
            r(:,i) = r(:,i-1) + longitud*orientacion
        end do
        ! Todos los N utilizados son pares.
        do i = 1,n,2
            do j = 1,3
                call normales(z1,z2)
                p(j,i) = sqrt(m*kBT)*z1
                p(j,i+1) = sqrt(m*kBT)*z2
            end do
        end do
    end subroutine inicializar

    subroutine paso_verlet(n,r,p,f,ue,f_anterior,ruido)
        integer, intent(in) :: n
        real(dp), intent(inout) :: r(3,n),p(3,n),f(3,n),ue
        real(dp), intent(out) :: f_anterior(3,n),ruido(3,n)
        integer :: i,j
        real(dp) :: z1,z2
        do i = 1,n,2
            do j = 1,3
                call normales(z1,z2)
                ruido(j,i) = amplitud_ruido*z1
                ruido(j,i+1) = amplitud_ruido*z2
            end do
        end do
        f_anterior = f
        r = r + b*h*p/m + b*h*h*f_anterior/(2.0_dp*m) + &
            b*h*ruido/(2.0_dp*m)
        call fuerzas(r,f,ue,l0,ke)
        p = a*p + h*(a*f_anterior+f)/2.0_dp + b*ruido
    end subroutine paso_verlet

    subroutine medir_geometria(r,rg2,ree2)
        real(dp), intent(in) :: r(:,:)
        real(dp), intent(out) :: rg2,ree2
        real(dp) :: rcm(3)
        integer :: n,i
        n = size(r,2)
        rcm = sum(r,dim=2)/real(n,dp)
        rg2 = 0.0_dp
        do i = 1,n
            rg2 = rg2 + sum((r(:,i)-rcm)**2)
        end do
        rg2 = rg2/real(n,dp)
        ree2 = sum((r(:,n)-r(:,1))**2)
    end subroutine medir_geometria

    subroutine simular_cadena(n)
        integer, intent(in) :: n
        real(dp), allocatable :: r(:,:),p(:,:),f(:,:),f_anterior(:,:),ruido(:,:)
        real(dp) :: ue,rg2,ree2,suma_rg2,suma_ree2
        real(dp) :: bloque_rg2,bloque_ree2,mitad1_rg2,mitad2_rg2
        integer :: paso,unidad,n_muestras,n_bloque,n_mitad1,n_mitad2
        character(len=240) :: nombre

        allocate(r(3,n),p(3,n),f(3,n),f_anterior(3,n),ruido(3,n))
        call inicializar(n,r,p)
        call fuerzas(r,f,ue,l0,ke)

        print *, 'N=',n,' empezando termalizacion'
        do paso = 1,pasos_termal
            call paso_verlet(n,r,p,f,ue,f_anterior,ruido)
        end do

        call ruta_salida(n,nombre)
        open(newunit=unidad,file=trim(nombre),status='new',action='write')
        write(unidad,'(A)') &
             '# t_medicion  Rg2_acumulado  Ree2_acumulado  Rg2_bloque  Ree2_bloque'
        suma_rg2=0.0_dp
        suma_ree2=0.0_dp
        bloque_rg2=0.0_dp
        bloque_ree2=0.0_dp
        mitad1_rg2=0.0_dp
        mitad2_rg2=0.0_dp
        n_muestras=0
        n_bloque=0
        n_mitad1=0
        n_mitad2=0

        print *, 'N=',n,' empezando medicion'
        do paso=1,pasos_medicion
            call paso_verlet(n,r,p,f,ue,f_anterior,ruido)
            if (mod(paso,cada_muestra)==0) then
                call medir_geometria(r,rg2,ree2)
                n_muestras=n_muestras+1
                n_bloque=n_bloque+1
                suma_rg2=suma_rg2+rg2
                suma_ree2=suma_ree2+ree2
                bloque_rg2=bloque_rg2+rg2
                bloque_ree2=bloque_ree2+ree2
                if (paso <= pasos_medicion/2) then
                    mitad1_rg2=mitad1_rg2+rg2
                    n_mitad1=n_mitad1+1
                else
                    mitad2_rg2=mitad2_rg2+rg2
                    n_mitad2=n_mitad2+1
                end if
            end if
            if (mod(paso,cada_guardado)==0) then
                write(unidad,'(5(ES21.12,1X))') &
                    real(paso,dp)*h, &
                    suma_rg2/real(n_muestras,dp), &
                    suma_ree2/real(n_muestras,dp), &
                    bloque_rg2/real(n_bloque,dp), &
                    bloque_ree2/real(n_bloque,dp)
                bloque_rg2=0.0_dp
                bloque_ree2=0.0_dp
                n_bloque=0
            end if
        end do
        close(unidad)
        print *, 'N=',n,' <Rg2>=',suma_rg2/real(n_muestras,dp), &
                 ' <Ree2>=',suma_ree2/real(n_muestras,dp)
        print *, '  <Rg2> primera mitad=',mitad1_rg2/real(n_mitad1,dp)
        print *, '  <Rg2> segunda mitad=',mitad2_rg2/real(n_mitad2,dp)
        print *, 'Guardado: ',trim(nombre)
        deallocate(r,p,f,f_anterior,ruido)
    end subroutine simular_cadena

end program main



