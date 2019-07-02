PROGRAM correlation
  IMPLICIT NONE

  ! Parameters/Variables
  ! Difference between two energy level \alpha = \omega_{gf} - \omega_{ge}
  REAL(KIND=8), PARAMETER :: alpha = -120
  ! Drive strength
  REAL(KIND=8), PARAMETER :: omega = 0.000001
  ! Relative driving strength
  REAL(KIND=8), PARAMETER :: xi = 1.0
  ! Decay Rate
  REAL(KIND=8), PARAMETER :: gamma = 1.0

  ! Detuning stuff
  ! Square root part
  REAL(KIND=8), PARAMETER :: DEL_sq = SQRT(((0.5 * alpha)**2) + 2*((0.5*Omega)**2)*((xi**2) - 1))
  ! Shifted single photon resonance
  REAL(KIND=8), PARAMETER :: DEL_1 = (-0.25 * alpha) + DEL_sq
  ! Shifted two-photon resonance
  REAL(KIND=8), PARAMETER :: DEL_2 = (-0.25 * alpha) - 0.5*DEL_sq
  ! Drive detuning
  REAL(KIND=8), PARAMETER :: delta =  0

  ! Operators
  ! Hamiltonian for the system (Wallraff group)
  REAL(KIND=8), DIMENSION(3,3) :: H
  ! Simplified Hamiltonian
  REAL(KIND=8), DIMENSION(3,3) :: H_s
  ! Energy level shifts
  REAL(KIND=8) :: d_1, d_2, d_g, d_e, d_f
  ! Effective drive strength
  REAL(KIND=8) :: o_eff
  ! Density Operators rho and steady state density operator
  COMPLEX(KIND=8), DIMENSION(3,3) :: rho, rho_ss
  ! Raising and Lowering Operators
  REAL(KIND=8), DIMENSION(3,3) :: sigmap, sigmam
  ! sigmap * sigmam
  REAL(KIND=8), DIMENSION(3,3) :: sigmapm
  ! sigmap * rho_ss
  COMPLEX(KIND=8), DIMENSION(3,3) :: sigrho
  ! sigrho copy for calculating TRACE(MATMUL(sigmap, sigrho))
  COMPLEX(KIND=8), DIMENSION(3,3) :: corr
  ! Density Operators for simplified version
  COMPLEX(KIND=8), DIMENSION(3,3) :: rho_s, rho_s_ss, sigrho_s, corr_simple

  ! Time step stuff
  ! Time step
  REAL(KIND=8), PARAMETER :: dt = 0.001
  ! Max time
  REAL(KIND=8), PARAMETER :: steady_time = 100
  ! Max number of time steps for tau calculations
  INTEGER, PARAMETER :: steady_steps = INT(steady_time / dt)
  ! Maximum tau time to calculate correlation function for
  REAL(KIND=8), PARAMETER :: tau_max = 50
  ! Max number of tau steps
  INTEGER, PARAMETER :: tau_steps = INT(tau_max / dt)

  ! Useful Thingies
  ! Useful integer place holders
  INTEGER :: l
  ! Complex i = SQRT(-1)
  COMPLEX(KIND=8), PARAMETER :: i = CMPLX(0,1)
  ! Runge-Kutta 4th Order Vectors/Matrices
  COMPLEX(KIND=8), DIMENSION(3,3) :: k1, k2, k3, k4
  ! Runge-Kutta vectors for the simple system
  COMPLEX(KIND=8), DIMENSION(3,3) :: ks1, ks2, ks3, ks4
  ! Temporaral value, temporal trace
  REAL(KIND=8) :: trace, ss_trace, trace_simple, ss_trace_simple
  ! Filenames
  CHARACTER(len=26) :: filename = "./data_files/corr_full.txt"
  CHARACTER(len=28) :: filename_simple = "./data_files/corr_simple.txt"

  ! Initialising matrices                 <g|<e|<f|
  ! In the atom basis the kets are: |g> = (1, 0, 0)^T
  !                                 |e> = (0, 1, 0)^T
  !                                 \f> = (0, 0, 1)^T

  ! Hamiltonian
  H = 0
  H(1,2) = 0.5 * omega
  H(2,1) = 0.5 * omega
  H(2,2) = -(0.5 * alpha + delta)
  H(2,3) = 0.5 * xi * omega
  H(3,2) = 0.5 * xi * omega
  H(3,3) = -2.0 * delta

  ! Simplified Hamiltonian
  IF (delta /= -0.5 * alpha) THEN
    o_eff = 2 * xi * ((0.5 * omega)**2) / (0.5*alpha + delta)
    d_1 = ((0.5 * omega) ** 2.0) / (0.5 * alpha + delta)
    d_2 = ((0.5 * xi * omega) ** 2.0) / (0.5 * alpha - delta)

    d_g = d_1
    d_f = d_2
    d_e = -1.0 * (d_1 + d_2)

    H_s = 0
    H_s(1,1) = d_g
    H_s(1,2) = 0.0
    H_s(1,3) = 0.5 * o_eff
    H_s(2,1) = 0.0
    H_s(2,2) = -(0.5 * alpha + delta) + d_e
    ! H_s(2,2) = 0.0
    H_s(2,3) = 0.0
    H_s(3,1) = 0.5 * o_eff
    H_s(3,2) = 0.0
    H_s(3,3) = (-2.0 * delta) + d_f
  ELSE
    H_S = 0
  END IF

  ! Raising and Lowering Operators
  sigmam = 0
  sigmam(1,2) = 1.0
  sigmam(2,3) = xi

  sigmap = 0
  sigmap(2,1) = sigmam(1,2)
  sigmap(3,2) = sigmam(2,3)

  ! S_+ x S_-
  sigmapm = 0
  sigmapm = MATMUL(sigmap, sigmam)

  ! Density rho initially in the ground state <g|rho|g> = 1
  rho = 0
  rho(1,1) = 1.0
  rho_s = 0
  rho_s(1,1) = 1.0

  ! Runge-Kutta 4-th Order Vectors
  k1 = 0
  k2 = 0
  k3 = 0
  k4 = 0

  ks1 = 0
  ks2 = 0
  ks3 = 0
  ks4 = 0

  ! Open the file used to write data
  OPEN(UNIT = 1, file = filename, STATUS = 'replace', ACTION = 'write')!
  OPEN(UNIT = 2, FILE = filename_simple, STATUS = 'replace', ACTION = 'write')

  ! Wrote the values of the parameters to the first 4 lines. The correlation
  !and times will start on the 6th line.
  WRITE(1,*), omega
  WRITE(1,*), delta
  WRITE(1,*), xi
  WRITE(1,*), alpha
  WRITE(1,*), ' '
  WRITE(2,*), omega
  WRITE(2,*), delta
  WRITE(2,*), xi
  WRITE(2,*), alpha
  WRITE(2,*), ' '

  ! Solve the master equation for a long time to find the system's steady state
  ! density matrix
  DO l=0,steady_steps   ! dt * steps = 0.01 * 1E5 = 1000 units
    ! Full model
    k1 = -i * dt * ( MATMUL(H, rho) - MATMUL(rho, H) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, rho), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, rho) - &
       & 0.5 * gamma * dt * MATMUL(rho, sigmapm)

    k2 = -i * dt * ( MATMUL(H, (rho + 0.5*k1)) - MATMUL((rho + 0.5*k1), H) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (rho + 0.5*k1)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (rho + 0.5*k1)) - &
       & 0.5 * gamma * dt * MATMUL((rho + 0.5*k1), sigmapm)

    k3 = -i * dt * ( MATMUL(H, (rho + 0.5*k2)) - MATMUL((rho + 0.5*k2), H) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (rho + 0.5*k2)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (rho + 0.5*k2)) - &
       & 0.5 * gamma * dt * MATMUL((rho + 0.5*k2), sigmapm)

    k4 = -i * dt * ( MATMUL(H, (rho + k3)) - MATMUL((rho + k3), H) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (rho + k3)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (rho + k3)) - &
       & 0.5 * gamma * dt * MATMUL((rho + k3), sigmapm)

    rho = rho + (1.0 / 6.0) * (k1 + 2.0*(k2 + k3) + k4)

    ! Simple model
    ks1 = -i * dt * ( MATMUL(H_s, rho_s) - MATMUL(rho_s, H_s) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, rho_s), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, rho_s) - &
       & 0.5 * gamma * dt * MATMUL(rho_s, sigmapm)

    ks2 = -i * dt * ( MATMUL(H_s, (rho_s + 0.5*k1)) - MATMUL((rho_s + 0.5*ks1), H_s) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (rho_s + 0.5*ks1)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (rho_s + 0.5*ks1)) - &
       & 0.5 * gamma * dt * MATMUL((rho_s + 0.5*ks1), sigmapm)

    ks3 = -i * dt * ( MATMUL(H_s, (rho_s + 0.5*ks2)) - MATMUL((rho_s + 0.5*ks2), H_s) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (rho_s + 0.5*ks2)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (rho_s + 0.5*ks2)) - &
       & 0.5 * gamma * dt * MATMUL((rho + 0.5*ks2), sigmapm)

    ks4 = -i * dt * ( MATMUL(H_s, (rho_s + ks3)) - MATMUL((rho_s + ks3), H_s) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (rho_s + ks3)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (rho_s + ks3)) - &
       & 0.5 * gamma * dt * MATMUL((rho_s + ks3), sigmapm)

    rho_s = rho_s + (1.0 / 6.0) * (ks1 + 2.0*(ks2 + ks3) + ks4)
  END DO

  rho_ss = rho
  rho_s_ss = rho_s

  ! Re-initialising the 4-th order vectors
  k1 = 0
  k2 = 0
  k3 = 0
  k4 = 0

  ks1 = 0
  ks2 = 0
  ks3 = 0
  ks4 = 0

  sigrho = 0
  sigrho = MATMUL(sigmapm, rho_ss)
  ss_trace = sigrho(1,1) + sigrho(2,2) + sigrho(3,3)
  !PRINT*, ss_trace

  ! Simplified model
  sigrho_s = 0
  sigrho_s = MATMUL(sigmapm, rho_s_ss)
  ss_trace_simple = sigrho_s(1,1) + sigrho_s(2,2) + sigrho_s(3,3)

  sigrho = 0
  ! Correlation functino given by Tr[ sigmap x exp(L*tau) x (sigmap x rho_ss) ]
  ! simgmam x rho_ss
  sigrho = MATMUL(MATMUL(sigmam, rho_ss), sigmap)

  ! Simplified model
  sigrho_s = 0
  sigrho_s = MATMUL(MATMUL(sigmam, rho_s_ss), sigmap)

  ! Calculate the first order correlation function lim t-> inft <sigmap(t+tau)sigmam(t)>
  DO l=0,tau_steps
    ! Instead of using sigrho to calculate the trace, create a copy so that
    ! sigrho remains untouched and is only used in the master equation evolution
    corr = MATMUL(sigmapm, sigrho)
    ! Calculate the trace
    trace = 0
    trace = (corr(1,1) + corr(2,2) + corr(3,3))
    trace = trace / (ss_trace ** 2)
    ! Write the results to file
    WRITE(1,*) l*dt, REAL(trace)

    ! Simplified model
    corr_simple = MATMUL(sigmapm, sigrho_s)
    trace_simple = 0
    trace_simple = (corr_simple(1,1) + corr_simple(2,2) +corr_simple(3,3))
    trace_simple = trace_simple / (ss_trace_simple ** 2)
    WRITE(2,*) l*dt, REAL(trace_simple)!, IMAG(trace_simple)

    ! Evolve the system up to a time tau
    k1 = -i * dt * ( MATMUL(H, sigrho) - MATMUL(sigrho, H) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, sigrho), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, sigrho) - &
       & 0.5 * gamma * dt * MATMUL(sigrho, sigmapm)

    k2 = -i * dt * ( MATMUL(H, (sigrho + 0.5*k1)) - MATMUL((sigrho + 0.5*k1), H) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (sigrho + 0.5*k1)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (sigrho + 0.5*k1)) - &
       & 0.5 * gamma * dt * MATMUL((sigrho + 0.5*k1), sigmapm)

    k3 = -i * dt * ( MATMUL(H, (sigrho + 0.5*k2)) - MATMUL((sigrho + 0.5*k2), H) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (sigrho + 0.5*k2)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (sigrho + 0.5*k2)) - &
       & 0.5 * gamma * dt * MATMUL((sigrho + 0.5*k2), sigmapm)

    k4 = -i * dt * ( MATMUL(H, (sigrho + k3)) - MATMUL((sigrho + k3), H) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (sigrho + k3)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (sigrho + k3)) - &
       & 0.5 * gamma * dt * MATMUL((sigrho + k3), sigmapm)

    sigrho = sigrho + (1.0 / 6.0) * (k1 + 2.0*(k2 + k3) + k4)

    ! Simple model
    ks1 = -i * dt * ( MATMUL(H_s, sigrho_s) - MATMUL(sigrho_s, H_s) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, sigrho_s), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, sigrho_s) - &
       & 0.5 * gamma * dt * MATMUL(sigrho_s, sigmapm)

    ks2 = -i * dt * ( MATMUL(H_s, (sigrho_s + 0.5*k1)) - MATMUL((sigrho_s + 0.5*ks1), H_s) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (sigrho_s + 0.5*ks1)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (sigrho_s + 0.5*ks1)) - &
       & 0.5 * gamma * dt * MATMUL((sigrho_s + 0.5*ks1), sigmapm)

    ks3 = -i * dt * ( MATMUL(H_s, (sigrho_s + 0.5*ks2)) - MATMUL((sigrho_s + 0.5*ks2), H_s) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (sigrho_s + 0.5*ks2)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (sigrho_s + 0.5*ks2)) - &
       & 0.5 * gamma * dt * MATMUL((sigrho_s + 0.5*ks2), sigmapm)

    ks4 = -i * dt * ( MATMUL(H_s, (sigrho_s + ks3)) - MATMUL((sigrho_s + ks3), H_s) ) + &
       & gamma * dt * MATMUL(MATMUL(sigmam, (sigrho_s + ks3)), sigmap) - &
       & 0.5 * gamma * dt * MATMUL(sigmapm, (sigrho_s + ks3)) - &
       & 0.5 * gamma * dt * MATMUL((sigrho_s + ks3), sigmapm)

    sigrho_s = sigrho_s + (1.0 / 6.0) * (ks1 + 2.0*(ks2 + ks3) + ks4)

  END DO

END PROGRAM correlation
