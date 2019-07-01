PROGRAM correlation
  IMPLICIT NONE

  ! Parameters/Variables
  REAL(KIND=8), PARAMETER :: alpha = -120
  ! Drive strength max
  REAL(KIND=8), PARAMETER :: omega =  40.0
  ! Drive detuning (calculate two-photon resonance) centre
  REAL(KIND=8), PARAMETER :: delta_centre  = 0.0
  ! Drive detuning displacement
  INTEGER, PARAMETER :: delta_displace = 20.0
  ! Relative driving strength
  REAL(KIND=8), PARAMETER :: xi = 1.5
  ! Decay Rate
  REAL(KIND=8), PARAMETER :: gamma = 1.0

  ! Operators
  ! Hamiltonian for the system
  COMPLEX(KIND=8), DIMENSION(3,3) :: H
  ! Density Operators rho and steady state density operator
  COMPLEX(KIND=8), DIMENSION(3,3) :: rho, rho_ss
  ! Raising and Lowering Operators
  COMPLEX(KIND=8), DIMENSION(3,3) :: sigmap, sigmam
  ! sigmap * sigmam
  COMPLEX(KIND=8), DIMENSION(3,3) :: sigmapm
  ! sigmap * rho_ss
  COMPLEX(KIND=8), DIMENSION(3,3) :: sigrho
  ! sigrho copy for calculating TRACE(MATMUL(sigmap, sigrho))
  COMPLEX(KIND=8), DIMENSION(3,3) :: corr

  ! Time thingies
  ! Time-step
  REAL, PARAMETER :: dt = 0.01
  ! Steady state time
  REAL, PARAMETER :: steady_time = 100
  ! Number of time steps for steady state
  INTEGER, PARAMETER :: steady_steps = INT(steady_time / dt)
  ! Max time
  REAL, PARAMETER :: tau_max = 20
  ! Max number of steps
  INTEGER, PARAMETER :: max_steps = INT(tau_max / dt)

  ! Parameter array
  ! Dimension of arrays
  INTEGER, PARAMETER :: N = 201
  ! Halfway integer
  INTEGER, PARAMETER :: Nhalf = INT(0.5 * (N-1))
  ! List of detunings delta
  REAl(KIND=8), DIMENSION(-Nhalf:Nhalf) :: delta_list
  ! Paramter to be called in loop
  REAL(KIND=8) :: delta

  ! Useful Thingies
  ! Useful integer place holders
  INTEGER :: j, k, l
  ! Complex i = SQRT(-1)
  COMPLEX(KIND=8), PARAMETER :: i = CMPLX(0,1)
  ! Runge-Kutta 4th Order Vectors/Matrices
  COMPLEX(KIND=8), DIMENSION(3,3) :: k1, k2, k3, k4
  ! Temporaral value, temporal trace
  COMPLEX(KIND=8) :: trace
  ! Steady state <\Sigma^{\dagger} \Sigma>
  REAL(KIND=8) ::  sigmapm_trace
  ! 1 / 6 as a constant
  REAL(KIND=8), PARAMETER :: xis = 1.0 / 6.0

  ! Filename stuff
  CHARACTER(len=31) :: filename_time = "./data_files/map_delta/time.txt"
  ! Filename for parameters
  CHARACTER(len=32) :: filename_parameters = "./data_files/map_delta/param.txt"
  ! Filename directory
  CHARACTER(len=27) :: file_directory = "./data_files/map_delta/data"
  ! Number of file for file name
  CHARACTER(len=4) :: file_number
  ! Filename extension
  CHARACTER(len=4) :: file_extension = ".txt"
  ! Combined filename
  CHARACTER(len=35) :: filename_total

  ! Delta array
  delta_list = 0
  DO j=-Nhalf,Nhalf
      delta_list(j) = delta_centre + (((1.0*j)/(1.0*Nhalf)) * delta_displace)
  END DO

  ! Initialising matrices
  ! Raising and Lowering
  sigmam = 0
  sigmam(1,2) = 1.0
  sigmam(2,3) = xi

  sigmap = 0
  sigmap(2,1) = sigmam(1,2)
  sigmap(3,2) = sigmam(2,3)

  sigmapm = 0
  sigmapm = MATMUL(sigmap, sigmam)

  ! Writing time values to file
  OPEN(UNIT=1, file=filename_time, STATUS='replace', ACTION='write')

  DO l=0,max_steps
    WRITE(1,*), l*dt
  END DO

  ! Writing parameter values to file
  OPEN(UNIT=3, file=filename_parameters, STATUS='replace', ACTION='write')
  WRITE(3,*), omega
  WRITE(3,*), xi
  WRITE(3,*), alpha
  WRITE(3,*), ' '

  DO j=-Nhalf,Nhalf
    ! Open a file for each value of j for the different deltas
    WRITE(file_number, '(I4.4)') j+Nhalf

    filename_total = file_directory // file_number // file_extension
    OPEN(UNIT = 2, file = filename_total, STATUS = 'replace', ACTION = 'write')

    delta = delta_list(j)
    WRITE(3,*), delta

    ! Hamiltonian
    H = 0
    H(1,2) = 0.5 * omega
    H(2,1) = 0.5 * omega
    H(2,2) = -(0.5 * alpha + delta)
    H(2,3) = 0.5 * xi * omega
    H(3,2) = 0.5 * xi * omega
    H(3,3) = -2.0 * delta

    ! Density rho initially in the ground state <g|rho|g> = 1 ?
    rho = 0
    rho(1,1) = 1.0

    ! Runge-Kutta
    k1 = 0
    k2 = 0
    k3 = 0
    k4 = 0

    ! Solve the master equation for a long time to find the system's steady state
    ! density matrix
    DO k=0,steady_steps
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

      rho = rho + xis * (k1 + 2.0*(k2 + k3) + k4)
    END DO

    rho_ss = 0
    rho_ss = rho

    ! Find steady state expectation <\Sigma^{\dagger} \Sigma>
    ! sigrho = MATMUL(sigmapm, rho_ss)
    sigrho = MATMUL(rho_ss, sigmapm)
    sigmapm_trace = 0
    DO k=1,3
      sigmapm_trace = sigmapm_trace + sigrho(k,k)
    END DO

    ! Reinitialising Runge-Kutta matrices
    k1 = 0
    k2 = 0
    k3 = 0
    k4 = 0

    ! Correlation functino given by Tr[ sigmap x exp(L*tau) x (sigmap x rho_ss) ]
    ! simgmam x rho_ss
    sigrho = MATMUL(rho_ss, sigmap)

    ! Calculate the first order correlation function lim t-> inft <sigmap(t+tau)sigmam(t)>
    DO k=0,max_steps
      ! Instead of using sigrho to calculate the trace, create a copy so that
      ! sigrho remains untouched and is only used in the master equation evolution
      corr = MATMUL(sigmam, sigrho)

      ! Calculate the trace
      trace = 0
      trace = corr(1,1) + corr(2,2) + corr(3,3)
      IF (sigmapm_trace /= 0) THEN
        trace = trace / REAL(sigmapm_trace)
      END IF

      ! Write the results to file
      WRITE(2,*) REAL(trace), IMAG(trace)

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

      sigrho = sigrho + xis * (k1 + 2.0*(k2 + k3) + k4)
    END DO

    PRINT*, 100.0 * (j + Nhalf) / (N), '% complete' !.Delta = ', delta

  END DO

END PROGRAM correlation
