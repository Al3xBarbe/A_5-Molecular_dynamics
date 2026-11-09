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

  use iso_fortran_env, only: real64, int64, output_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: dp = real64

  ! ==========================================================
  ! 1. PARAMETROS FISICOS
  ! ==========================================================

  integer, parameter :: sizes(10) = &
       [4,8,12,16,20,24,32,40,48,64]

  real(dp), parameter :: dt = 0.001_dp
  real(dp), parameter :: ke = 1000.0_dp
  real(dp), parameter :: kbt = 1.0_dp
  real(dp), parameter :: l0 = 1.0_dp
  real(dp), parameter :: masa = 1.0_dp
  real(dp), parameter :: eta = 1.0_dp

  integer, parameter :: nbloques = 20

  ! Una muestra cada 50 unidades de tiempo
  integer(int64), parameter :: salto_normal = 50000_int64

  character(len=*), parameter :: ruta_normal = &
       'Py/plots/OBJETIVO_3/Datos/k1000/EQUIP_E'

  ! ==========================================================
  ! 2. VARIABLES
  ! ==========================================================

  integer :: n, ix, j, ios, nsem
  integer :: argumento_n
  integer :: udatos, ubloques, umarca
  integer, allocatable :: semillas(:)

  integer(int64) :: nterm, nmed, salto
  integer(int64) :: paso, paso_med
  integer(int64) :: nmedidas, por_bloque
  integer(int64) :: cont_bloque, num_bloque

  real(dp) :: tterm, tmed
  real(dp) :: a, b, ruido
  real(dp) :: u, kin, rg, ree, rg2, ree2

  real(dp) :: sk, sv, srg, sree, srg2, sree2
  real(dp) :: bk, bv, brg, bree, brg2, bree2

  real(dp), allocatable :: r(:,:), rnew(:,:)
  real(dp), allocatable :: v(:,:), vnew(:,:)
  real(dp), allocatable :: f(:,:), fnew(:,:), beta(:,:)

  character(len=250) :: carpeta
  character(len=250) :: archivo_d, archivo_b, archivo_m
  character(len=80) :: arg
  character(len=20) :: ns

  logical :: existe, prueba
  real(dp) :: spare
  logical :: tiene_spare

  ! ==========================================================
  ! 3. ARGUMENTOS DE EJECUCION
  ! ==========================================================

  ! Sin argumentos: simula los 10 tamaños.
  ! Con un N: simula solamente ese tamaño.
  ! Con "prueba": realiza una prueba corta.

  argumento_n = 0
  prueba = .false.

  if (command_argument_count() >= 1) then

     call get_command_argument(1, arg)

     if (trim(arg) == 'prueba') then
        prueba = .true.
     else
        read(arg, *, iostat=ios) argumento_n

        if (ios /= 0) error stop 'Argumento no valido'

        if (.not. any(sizes == argumento_n)) then
           error stop 'N no valido'
        end if
     end if

  end if

  if (command_argument_count() >= 2) then

     call get_command_argument(2, arg)

     if (trim(arg) /= 'prueba') then
        error stop 'Segundo argumento debe ser prueba'
     end if

     prueba = .true.

  end if

  if (prueba) then
     carpeta = ruta_normal//'/PRUEBA'
  else
     carpeta = ruta_normal
  end if

  call execute_command_line( &
       'mkdir -p "'//trim(carpeta)//'"', exitstat=ios)

  if (ios /= 0) error stop 'No se puede crear carpeta'

  ! ==========================================================
  ! 4. INTEGRADOR GJF
  ! ==========================================================

  a = (1.0_dp - eta*dt/(2.0_dp*masa)) / &
      (1.0_dp + eta*dt/(2.0_dp*masa))

  b = 1.0_dp / &
      (1.0_dp + eta*dt/(2.0_dp*masa))

  ruido = sqrt(2.0_dp*eta*kbt*dt)

  call random_seed(size=nsem)
  allocate(semillas(nsem))

  print *, '============================================'
  print *, 'OBJETIVO 3 - SIMULACIONES ke = 1000'
  print *, 'h = ', dt
  print *, 'N = 4,8,12,16,20,24,32,40,48,64'
  print *, '============================================'

  flush(output_unit)

  ! ==========================================================
  ! 5. BUCLE SOBRE TAMAÑOS DE CADENA
  ! ==========================================================

  do ix = 1, size(sizes)

     n = sizes(ix)

     if (argumento_n /= 0 .and. n /= argumento_n) cycle

     ! -------------------------------------------------------
     ! TIEMPOS ACORDADOS
     ! -------------------------------------------------------

     tterm = 1000000.0_dp

     select case (n)

     case (4,8,12,16,20)
        tmed = 3000000.0_dp

     case (24)
        tmed = 2750000.0_dp

     case (32)
        tmed = 2500000.0_dp

     case (40)
        tmed = 2250000.0_dp

     case (48,64)
        tmed = 2000000.0_dp

     end select

     salto = salto_normal

     ! Prueba corta con los mismos parametros fisicos
     if (prueba) then
        tterm = 100.0_dp
        tmed = 200.0_dp
        salto = 1000_int64
     end if

     nterm = nint(tterm/dt, kind=int64)
     nmed = nint(tmed/dt, kind=int64)

     if (mod(nterm, salto) /= 0_int64) then
        error stop 'Termalizacion no divisible por salto'
     end if

     if (mod(nmed, salto*nbloques) /= 0_int64) then
        error stop 'Medicion no divisible en bloques'
     end if

     por_bloque = nmed / (salto*nbloques)

     ! -------------------------------------------------------
     ! NOMBRES DE ARCHIVOS
     ! -------------------------------------------------------

     write(ns, '(I0)') n

     archivo_d = trim(carpeta)// &
       '/Energia_Rg_Ree_N'//trim(ns)//'_ke1000_ULTRALARGA.txt'

     archivo_b = trim(carpeta)// &
       '/Bloques_N'//trim(ns)//'_ke1000_ULTRALARGA.txt'

     archivo_m = trim(carpeta)// &
       '/COMPLETADO_N'//trim(ns)//'_ke1000_ULTRALARGA.flag'

     ! -------------------------------------------------------
     ! PROTEGER LOS ARCHIVOS YA FINALIZADOS
     ! -------------------------------------------------------

     inquire(file=trim(archivo_m), exist=existe)

     if (existe) then
        print '(A,I0,A)', &
          'N=', n, ': ya completado, se omite.'
        cycle
     end if

     inquire(file=trim(archivo_d), exist=existe)

     if (existe) then
        error stop 'Existe un TXT parcial: revisarlo'
     end if

     inquire(file=trim(archivo_b), exist=existe)

     if (existe) then
        error stop 'Existe archivo parcial de bloques'
     end if

     ! -------------------------------------------------------
     ! SEMILLA REPRODUCIBLE POR CADENA
     ! -------------------------------------------------------

     do j = 1, nsem
        semillas(j) = 94717 + 1301*n + 41*j
     end do

     call random_seed(put=semillas)
     tiene_spare = .false.

     allocate(r(3,n), rnew(3,n))
     allocate(v(3,n), vnew(3,n))
     allocate(f(3,n), fnew(3,n), beta(3,n))

     ! Inicializacion de equilibrio
     call inicializar(n, r, v)
     call fuerzas(n, r, f, u)

     ! -------------------------------------------------------
     ! ABRIR ARCHIVOS
     ! -------------------------------------------------------

     open(newunit=udatos, &
       file=trim(archivo_d), &
       status='new', action='write', iostat=ios)

     if (ios /= 0) error stop 'Error abriendo datos'

     open(newunit=ubloques, &
       file=trim(archivo_b), &
       status='new', action='write', iostat=ios)

     if (ios /= 0) error stop 'Error abriendo bloques'

     write(udatos, '(A)') &
       '# t_med K V Rg Ree Rg2 Ree2 ' // &
       'K_acum V_acum Rg_acum Ree_acum Rg2_acum Ree2_acum'

     write(ubloques, '(A)') &
       '# bloque t_inicio t_fin K_b V_b ' // &
       'Rg_b Ree_b Rg2_b Ree2_b'

     ! -------------------------------------------------------
     ! INICIALIZAR PROMEDIOS
     ! -------------------------------------------------------

     sk=0.0_dp
     sv=0.0_dp
     srg=0.0_dp
     sree=0.0_dp
     srg2=0.0_dp
     sree2=0.0_dp

     bk=0.0_dp
     bv=0.0_dp
     brg=0.0_dp
     bree=0.0_dp
     brg2=0.0_dp
     bree2=0.0_dp

     nmedidas=0_int64
     cont_bloque=0_int64
     num_bloque=0_int64

     print *
     print '(A,I0,A,F12.0,A,F12.0)', &
       'INICIO N=', n, &
       ' / termalizacion=', tterm, &
       ' / medicion=', tmed

     flush(output_unit)

     ! =======================================================
     ! 6. INTEGRACION DE LANGEVIN
     ! =======================================================

     do paso = 1_int64, nterm+nmed

        ! Generar impulsos aleatorios
        call ruido_gauss(n, beta, ruido)

        ! Nueva posicion mediante GJF
        rnew = r + b*dt*v + &
          (b*dt*dt/(2.0_dp*masa))*f + &
          (b*dt/(2.0_dp*masa))*beta

        ! Fuerzas y energia potencial nuevas
        call fuerzas(n, rnew, fnew, u)

        ! Nueva velocidad interna GJF
        vnew = a*v + &
          (dt/(2.0_dp*masa))*(a*f+fnew) + &
          (b/masa)*beta

        ! ====================================================
        ! FASE DE MEDICION
        ! ====================================================

        if (paso > nterm) then

           paso_med = paso-nterm

           if (mod(paso_med, salto) == 0_int64) then

              ! Energia cinetica estimada con 2GJ
              kin = masa*sum((rnew-r)**2) / &
                    (2.0_dp*b*dt*dt)

              ! Radios y sus cuadrados
              call radios(n, rnew, rg, ree, rg2, ree2)

              ! Comprobacion de resultados no finitos
              if (.not. ieee_is_finite(kin) .or. &
                  .not. ieee_is_finite(u) .or. &
                  .not. ieee_is_finite(rg2) .or. &
                  .not. ieee_is_finite(ree2)) then
                 error stop 'Valor no finito en medicion'
              end if

              nmedidas = nmedidas+1_int64
              cont_bloque = cont_bloque+1_int64

              ! Sumas acumuladas
              sk = sk+kin
              sv = sv+u

              srg = srg+rg
              sree = sree+ree

              srg2 = srg2+rg2
              sree2 = sree2+ree2

              ! Sumas del bloque actual
              bk = bk+kin
              bv = bv+u

              brg = brg+rg
              bree = bree+ree

              brg2 = brg2+rg2
              bree2 = bree2+ree2

              ! ----------------------------------------------
              ! DATOS Y MEDIAS ACUMULADAS
              ! ----------------------------------------------

              write(udatos, '(13(ES23.14,1X))') &
                real(paso_med,dp)*dt, &
                kin, u, rg, ree, rg2, ree2, &
                sk/real(nmedidas,dp), &
                sv/real(nmedidas,dp), &
                srg/real(nmedidas,dp), &
                sree/real(nmedidas,dp), &
                srg2/real(nmedidas,dp), &
                sree2/real(nmedidas,dp)

              ! ----------------------------------------------
              ! MEDIAS DE 20 BLOQUES
              ! ----------------------------------------------

              if (cont_bloque == por_bloque) then

                 num_bloque = num_bloque+1_int64

                 write(ubloques, &
                   '(I5,1X,8(ES23.14,1X))') &
                   int(num_bloque), &
                   real(num_bloque-1_int64,dp)* &
                     tmed/real(nbloques,dp), &
                   real(num_bloque,dp)* &
                     tmed/real(nbloques,dp), &
                   bk/real(cont_bloque,dp), &
                   bv/real(cont_bloque,dp), &
                   brg/real(cont_bloque,dp), &
                   bree/real(cont_bloque,dp), &
                   brg2/real(cont_bloque,dp), &
                   bree2/real(cont_bloque,dp)

                 bk=0.0_dp
                 bv=0.0_dp
                 brg=0.0_dp
                 bree=0.0_dp
                 brg2=0.0_dp
                 bree2=0.0_dp

                 cont_bloque=0_int64

              end if

           end if

           ! Mostrar progreso de medicion
           if (mod(paso_med, nmed/10_int64) == 0_int64) then

              print '(A,I0,A,I0,A)', &
                'N=', n, ' / medicion ', &
                int(100_int64*paso_med/nmed), '%'

              flush(output_unit)

           end if

        ! ====================================================
        ! FASE DE TERMALIZACION
        ! ====================================================

        else

           if (mod(paso, nterm/4_int64) == 0_int64) then

              print '(A,I0,A,I0,A)', &
                'N=', n, ' / termalizacion ', &
                int(100_int64*paso/nterm), '%'

              flush(output_unit)

           end if

        end if

        ! Actualizar el estado del sistema
        r = rnew
        v = vnew
        f = fnew

     end do

     ! =======================================================
     ! 7. FINALIZACION
     ! =======================================================

     close(udatos, iostat=ios)
     if (ios /= 0) error stop 'Error cerrando datos'

     close(ubloques, iostat=ios)
     if (ios /= 0) error stop 'Error cerrando bloques'

     ! Marca de completado: solo al terminar correctamente
     open(newunit=umarca, &
       file=trim(archivo_m), &
       status='new', action='write', iostat=ios)

     if (ios /= 0) error stop 'No se crea marca'

     write(umarca, '(A,I0)') 'N=', n
     write(umarca, '(A,I0)') 'muestras=', nmedidas
     write(umarca, '(A,F14.3)') 'k_e=', ke
     write(umarca, '(A,F14.6)') 'h=', dt
     write(umarca, '(A,F14.0)') 'termalizacion=', tterm
     write(umarca, '(A,F14.0)') 'medicion=', tmed

     close(umarca)

     print *
     print '(A,I0,A,F14.6,A,F14.6)', &
       'FIN N=', n, &
       ' <K>=', sk/real(nmedidas,dp), &
       ' <V>=', sv/real(nmedidas,dp)

     print '(A,F14.6,A,F14.6)', &
       ' <Rg>=', srg/real(nmedidas,dp), &
       ' <Ree>=', sree/real(nmedidas,dp)

     flush(output_unit)

     deallocate(r, rnew, v, vnew)
     deallocate(f, fnew, beta)

  end do

  deallocate(semillas)

  print *
  print *, 'FIN DE LA TANDA ke = 1000'

contains

  ! ==========================================================
  ! 8. GENERADOR GAUSSIANO
  ! ==========================================================

  real(dp) function gauss() result(z)

    real(dp) :: x, y, s, factor

    if (tiene_spare) then

       z = spare
       tiene_spare = .false.

    else

       do

          call random_number(x)
          call random_number(y)

          x = 2.0_dp*x-1.0_dp
          y = 2.0_dp*y-1.0_dp

          s = x*x+y*y

          if (s > 0.0_dp .and. s < 1.0_dp) exit

       end do

       factor = sqrt(-2.0_dp*log(s)/s)

       z = x*factor
       spare = y*factor

       tiene_spare = .true.

    end if

  end function gauss

  ! ==========================================================
  ! 9. LONGITUD INICIAL CON JACOBIANO RADIAL
  ! ==========================================================

  real(dp) function longitud() result(l)

    real(dp) :: s2, xstar, cota, x, p, uale

    s2 = kbt/ke

    xstar = (-l0+sqrt(l0*l0+16.0_dp*s2))/2.0_dp

    cota = (l0+xstar)**2 * &
           exp(-xstar*xstar/(4.0_dp*s2))

    do

       x = sqrt(2.0_dp*s2)*gauss()
       l = l0+x

       if (l <= 0.0_dp) cycle

       p = l*l*exp(-x*x/(4.0_dp*s2))/cota

       call random_number(uale)

       if (uale <= min(1.0_dp,p)) exit

    end do

  end function longitud

  ! ==========================================================
  ! 10. INICIALIZACION TERMICA
  ! ==========================================================

  subroutine inicializar(nloc, pr, vr)

    integer, intent(in) :: nloc

    real(dp), intent(out) :: pr(3,nloc)
    real(dp), intent(out) :: vr(3,nloc)

    real(dp) :: z, phi, sint, l, ale
    integer :: i, jj

    pr(:,1) = 0.0_dp

    do jj = 2, nloc

       l = longitud()

       call random_number(ale)
       z = 2.0_dp*ale-1.0_dp

       call random_number(ale)
       phi = 2.0_dp*acos(-1.0_dp)*ale

       sint = sqrt(max(0.0_dp,1.0_dp-z*z))

       pr(:,jj) = pr(:,jj-1) + &
          l*[sint*cos(phi), sint*sin(phi), z]

    end do

    ! Maxwell-Boltzmann
    do jj = 1, nloc
       do i = 1, 3
          vr(i,jj) = sqrt(kbt/masa)*gauss()
       end do
    end do

  end subroutine inicializar

  ! ==========================================================
  ! 11. RUIDO ESTOCASTICO
  ! ==========================================================

  subroutine ruido_gauss(nloc, bt, amp)

    integer, intent(in) :: nloc

    real(dp), intent(in) :: amp
    real(dp), intent(out) :: bt(3,nloc)

    integer :: i, jj

    do jj = 1, nloc
       do i = 1, 3
          bt(i,jj) = amp*gauss()
       end do
    end do

  end subroutine ruido_gauss

  ! ==========================================================
  ! 12. FUERZAS Y ENERGIA POTENCIAL
  ! ==========================================================

  subroutine fuerzas(nloc, pr, ff, ep)

    integer, intent(in) :: nloc

    real(dp), intent(in) :: pr(3,nloc)
    real(dp), intent(out) :: ff(3,nloc), ep

    real(dp) :: d(3), l, c
    integer :: jj

    ff = 0.0_dp
    ep = 0.0_dp

    do jj = 1, nloc-1

       d = pr(:,jj+1)-pr(:,jj)

       l = sqrt(sum(d*d))

       if (l < 1.0e-12_dp) then
          error stop 'Enlace de longitud nula'
       end if

       c = ke*(l-l0)/l

       ff(:,jj) = ff(:,jj)+c*d
       ff(:,jj+1) = ff(:,jj+1)-c*d

       ep = ep + 0.5_dp*ke*(l-l0)**2

    end do

  end subroutine fuerzas

  ! ==========================================================
  ! 13. RADIO DE GIRO Y DISTANCIA ENTRE EXTREMOS
  ! ==========================================================

  subroutine radios(nloc, pr, rg, ree, rg2, ree2)

    integer, intent(in) :: nloc

    real(dp), intent(in) :: pr(3,nloc)

    real(dp), intent(out) :: rg, ree, rg2, ree2

    real(dp) :: cm(3)
    integer :: jj

    cm = sum(pr, dim=2)/real(nloc,dp)

    rg2 = 0.0_dp

    do jj = 1, nloc

       rg2 = rg2+sum((pr(:,jj)-cm)**2)

    end do

    rg2 = rg2/real(nloc,dp)

    ree2 = sum((pr(:,nloc)-pr(:,1))**2)

    rg = sqrt(rg2)
    ree = sqrt(ree2)

  end subroutine radios

end program main
