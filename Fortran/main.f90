!========================== generar_normal(s, a0, c0, m0, no1, no2, In2) ===========================================
!========================== generar_normal_txt(pasos, semilla, a, c, m, nombre_archivo) ===========================================
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


    subroutine generar_normal_txt(pasos, semilla, a, c, m, nombre_archivo)

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

    end subroutine generar_normal_txt

end module generador_normal

!========================== euler_maruyama(pasos,xo,po,m,k,beta_inv,nu) ===========================================
!========================== Runge_Kutta_2 (pasos,h,xo,po,m,k,beta_inv,nu) ===========================================
!========================== verlet_exp_est(pasos,h,xo,po,m,k,beta_inv,nu) ===========================================
module algoritmos_estocastios
    use generador_normal
    contains

    subroutine euler_maruyama(pasos,h,xo,po,m,k,beta_inv,nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        integer(kind=8) :: semilla, a, c, modulo
        real :: x, p, n1, n2, aux2, Ki, V, E
        integer(kind=8) :: aux
        integer :: i

        ! Parametros del generador
        semilla = 1256789
        a       = 164525
        c       = 1004223
        modulo  = 4294296

        x = xo
        p = po

        open(unit=3, file="euler_maruyama.txt", status="replace")
        Ki=0.5*p*p/m
        V=0.5*k*x*x
        E=Ki+V
        write(3,*) 0.0, x, p, Ki, V, E
        do i = 1, pasos/2

            call generar_normal(semilla, a, c, modulo, n1, n2, aux)
            semilla = aux

            ! Actualizacion de las variables usando Euler-Maruyama
            aux2 = x
            x = x + (p/m)*h 
            p = p - k*aux2*h - nu*p*h + sqrt(2.0*nu*m*beta_inv*h)*n1

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(3,*) (2*i-1)*h,x, p, Ki, V, E

            aux2 = x
            x = x + (p/m)*h 
            p = p - k*aux2*h - nu*p*h + sqrt(2.0*nu*m*beta_inv*h)*n2

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(3,*) (2*i)*h, x, p, Ki, V, E

        end do
        close(3)

    end subroutine euler_maruyama

    subroutine Runge_Kutta_2 (pasos, h, xo, po, m, k, beta_inv, nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        integer(kind=8) :: semilla, a, c, modulo
        real :: x, p, n1, n2, fx1, fx2, fg1, fg2, z_det, Ki, V, E
        integer(kind=8) :: aux
        integer :: i

        ! Parametros del generador
        semilla = 1256789
        a       = 164525
        c       = 1004223
        modulo  = 4294296

        x = xo
        p = po

        open(unit=4, file="runge_kutta_2.txt", status="replace")
        Ki=0.5*p*p/m
        V=0.5*k*x*x
        E=Ki+V
        write(4,*) 0.0,x, p, Ki, V, E
        z_det=sqrt(2.0*nu*m*beta_inv*h)
        do i = 1, pasos/2

            call generar_normal(semilla, a, c, modulo, n1, n2, aux)
            semilla = aux

            fx1 = (p+z_det*n1)/m
            fg1 = -k*x - nu*(p+z_det*n1)

            fx2 = (p + h*fg1)/2
            fg2 = -k*(x + h*fx1) - nu*(p + h*fg1)

            x = x + 0.5*h*(fx1 + fx2)
            p = p + 0.5*h*(fg1 + fg2) + z_det*n1

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(4,*) (2*i-1)*h, x, p, Ki, V, E

            fx1 = (p+z_det*n2)/m
            fg1 = -k*x - nu*(p+z_det*n2)

            fx2 = (p + h*fg1)/2
            fg2 = -k*(x + h*fx1) - nu*(p + h*fg1)

            x = x + 0.5*h*(fx1 + fx2)
            p = p + 0.5*h*(fg1 + fg2) + z_det*n2

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(4,*) (2*i)*h, x, p, Ki, V, E

        end do

        close(4)

    end subroutine Runge_Kutta_2

    subroutine verlet_exp_est(pasos, h, xo, po, m, k, beta_inv, nu)
        implicit none
        integer, intent(in) :: pasos
        real, intent(in) :: h, xo, po, m, k, beta_inv, nu
        integer(kind=8) :: semilla, a, c, modulo
        real :: x, p, n1, n2, z_det, verlet_a, verlet_b, aux2, Ki, V, E
        integer(kind=8) :: aux
        integer :: i

        ! Parametros del generador
        semilla = 1256789
        a       = 164525
        c       = 1004223
        modulo  = 4294296

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

            call generar_normal(semilla, a, c, modulo, n1, n2, aux)
            semilla = aux


            aux2=x
            x = x + verlet_b*h*p + (-verlet_b*h*h*x)/(2*m) + verlet_b*h*z_det*n1/(2*m)
            p = verlet_a*p - h*(verlet_a*k*aux2 + k*x)/(2*m) + verlet_b*z_det*n1/(2*m)

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(5,*) (2*i-1)*h, x, p, Ki, V, E

            aux2=x
            x = x + verlet_b*h*p + (-verlet_b*h*h*x)/(2*m) + verlet_b*h*z_det*n2/(2*m)
            p = verlet_a*p - h*(verlet_a*k*aux2 + k*x)/(2*m) + verlet_b*z_det*n2/(2*m)

            Ki=0.5*p*p/m
            V=0.5*k*x*x
            E=Ki+V
            write(5,*) (2*i)*h, x, p, Ki, V, E

        end do
        close(5)

    end subroutine verlet_exp_est

end module algoritmos_estocastios

program main

    use generador_normal
    use algoritmos_estocastios

    implicit none

    integer(kind=8) :: semilla, a, c, modulo
    integer :: pasos
    real :: h, x0, p0, m, k, beta_inv, nu

    ! Parametros del generador
    semilla = 1256789
    a       = 164525
    c       = 1004223
    modulo  = 4294296

    ! Parametros de la simulacion
    pasos = 100000
    h=0.01

    x0=5.0
    p0=0.0
    m=1.0
    k=1.0
    beta_inv=1.0
    nu=0.01 

    call euler_maruyama(pasos, h,x0,p0,m,k,beta_inv,nu)
    call Runge_Kutta_2(pasos, h, x0, p0, m, k, beta_inv, nu)
    call verlet_exp_est(pasos, h, x0, p0, m, k, beta_inv, nu)

    read(*,*)


end program main