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
        real :: x, p, n1, n2, z_det, verlet_a, verlet_b, aux
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


        open(unit=176, file="verlet_doblepozo.txt", status="replace")
        open(unit=177, file="estancia_pozo1.txt", status="replace")
        open(unit=178, file="estancia_pozo2.txt", status="replace")

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
                if (aux>0) then
                    write(178,*) Contador2*h
                else
                    write(177,*) Contador2*h
                end if
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
                if (aux>0.0) then
                    write(178,*) Contador2*h
                else
                    write(177,*) Contador2*h
                end if
                Contador2 = 0
            end if

        end do

        close(176)
        close(177)
        close(178)

    end subroutine verlet_doblepozo

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
        real(dp) :: rcm(3), rg2, rg2_m, ue
        integer :: i,j,ka

        open(unit=181, file="Radio_giro.txt", status="replace")

        z_det=sqrt(2.0_dp*eta*m*beta_inv*h)
        verlet_a = (1.0_dp-(eta*h)/2.0_dp)/(1.0_dp+(eta*h)/2.0_dp)
        verlet_b = 1.0_dp/(1.0_dp+(eta*h)/2.0_dp)
        call fuerzas(r,f,ue,l0,k)

        rcm = sum(r, dim=2) / size(r,2)
        rg2 = 0.0_dp
        do i = 1, size(r,2)
            rg2 = rg2 + sum((r(:,i) - rcm)**2)
        end do
        rg2 = rg2 / size(r,2)

        rg2_m=rg2

        write(181,*) 0.0_dp, rg2_m


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

            rg2_m=rg2_m + rg2

            if (mod(i,50) == 0) write(181,*) (2*i-1)*h, rg2_m/(2*i)

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

            rg2_m=rg2_m + rg2

            if (mod(i,50) == 0) write(181,*) 2*i*h, rg2_m/(2*i+1)

        end do

        close(181)

    end subroutine radio_giro

end module objetivo_3

program main

    use generador_normal
    use algoritmos_estocastios
    use objetivo_2
    use objetivo_3
    implicit none

!OBJETIVO 2

    integer :: tiempo, pasos
    real :: F,x0,p0,A, h, m, beta_inv, eta

    tiempo=500000
    h=0.001
    pasos=int(real(tiempo)/h)
    
    x0=1.0
    p0=0.0
    F=-0.80
    A=0.5
    m = 1.0
    beta_inv = 0.2
    eta = 0.4

!OBJETIVO 3

    !pasos=int(real(tiempo,dp)/h)
    !h=0.001_dp
    !integer :: i
    !integer, parameter :: n=4
    !real(dp) :: h, m, ke, beta_inv, eta, l0
    !real(dp) :: r(3,n), p(3,n)

    !tiempo=500000
    !l0=1.0_dp
    !do i=1, n
    !    r(:,i)=[real(i-1,dp)*l0, 0.0_dp, 0.0_dp]
    !end do
    !p=0.0_dp
    !m=1.0_dp
    !ke=100.0_dp
    !beta_inv=1.0_dp
    !eta=1.0_dp


    !call radio_giro(pasos,h,r,p,m,ke,beta_inv,eta,l0)
    call verlet_dp_fuerza(pasos,h,x0,p0,m,A,beta_inv,eta,F)

    print *, "Pulsa ENTER para salir..."
    read(*,*)
end program main