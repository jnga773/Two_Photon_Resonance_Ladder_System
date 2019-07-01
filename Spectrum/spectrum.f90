PROGRAM first_order_correlation
  IMPLICIT NONE

  ! Parameters/Variables
  REAL(KIND=8), PARAMETER :: alpha = -120.0
  ! Drive strength
  REAL(KIND=8), PARAMETER :: omega = 5.0
  ! Relative driving strength
  REAL(KIND=8), PARAMETER :: xi = 0.0
  ! Decay Rate
  REAL(KIND=8), PARAMETER :: gamma = 1.0

  ! Detuning stuff
  ! Square root part
  REAL(KIND=8), PARAMETER :: DEL_sq = SQRT(((0.5 * alpha)**2) + 2*((0.5*Omega)**2)*((xi**2) - 1))
  ! Shifted single photon resonance
  REAL(KIND=8), PARAMETER :: DEL_1 = (-0.25 * alpha) + 0.5 * DEL_sq
  ! Shifted two-photon resonance
  REAL(KIND=8), PARAMETER :: DEL_2 = (-0.25 * alpha) - 0.5 * DEL_sq
  ! Drive detuning
  REAL(KIND=8), PARAMETER :: delta = 60.0

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
  REAL, PARAMETER :: steady_time = 50
  ! Steady state steps
  INTEGER, PARAMETER :: steady_steps = INT(steady_time / dt)
  ! Max time
  REAL, PARAMETER :: t_max = 1000
  ! Max number of steps
  INTEGER, PARAMETER :: max_steps = INT(t_max / dt)

  ! Useful Thingies
  ! Useful integer place holders
  INTEGER :: k
  ! Complex i = SQRT(-1)
  COMPLEX(KIND=8), PARAMETER :: i = CMPLX(0,1)
  ! Runge-Kutta 4th Order Vectors/Matrices
  COMPLEX(KIND=8), DIMENSION(3,3) :: k1, k2, k3, k4
  ! Temporaral value, temporal trace
  COMPLEX(KIND=8) :: trace
  ! Steady state <\Sigma^{\dagger} \Sigma>
  COMPLEX(KIND=8) ::  sigmapm_trace, sigmam_trace, sigmap_trace
  ! Coherent and incoherent intensity ratio
  REAL(KIND=8) :: intensity_ratio
  ! 1 / 6 as a constant
  REAL(KIND=8), PARAMETER :: xis = 1.0 / 6.0
  !Filename for spectrum correlation
  CHARACTER(len=25) :: filename_spectrum = "./data_files/spectrum.txt"

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

  ! Runge-Kutta 4-th Order Vectors
  k1 = 0
  k2 = 0
  k3 = 0
  k4 = 0

  ! Open the file used to write data
  OPEN(UNIT = 1, file = filename_spectrum, STATUS = 'replace', ACTION = 'write')!

  ! Wrote the values of the parameters to the first 4 lines. The correlation
  !and times will start on the 6th line.
  WRITE(1,*), omega
  WRITE(1,*), delta
  WRITE(1,*), xi
  WRITE(1,*), alpha
  WRITE(1,*), ' '

  PRINT*, "Computing steady state"

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

  rho_ss = rho

  PRINT*, "Steady state found"

  ! Re-initialising the 4-th order vectors
  k1 = 0
  k2 = 0
  k3 = 0
  k4 = 0

  ! Find steady state expectation <\Sigma^{\dagger}>
  ! sigrho = MATMUL(sigmapm, rho_ss)
  sigrho = MATMUL(rho_ss, sigmap)
  sigmap_trace = 0
  DO k=1,3
    sigmap_trace = sigmap_trace + sigrho(k,k)
  END DO
  PRINT*, "<\Sigma^{\dagger}>_{ss} = ", sigmap_trace

  ! Find steady state expectation <\Sigma>
  ! sigrho = MATMUL(sigmapm, rho_ss)
  sigrho = MATMUL(rho_ss, sigmam)
  sigmam_trace = 0
  DO k=1,3
    sigmam_trace = sigmam_trace + sigrho(k,k)
  END DO
  PRINT*, "<\Sigma>_{ss} = ", sigmam_trace

  ! Find steady state expectation <\Sigma^{\dagger} \Sigma>
  ! sigrho = MATMUL(sigmapm, rho_ss)
  sigrho = MATMUL(rho_ss, sigmapm)
  sigmapm_trace = 0
  DO k=1,3
    sigmapm_trace = sigmapm_trace + sigrho(k,k)
  END DO
  PRINT*, "<\Sigma^{\dagger} \Sigma>_{ss} = ", sigmapm_trace

  ! Intensity ratio
  intensity_ratio = (REAL(sigmapm_trace) - REAL(sigmap_trace * sigmam_trace))
  intensity_ratio = intensity_ratio / REAL(sigmap_trace * sigmam_trace)
  PRINT*, "I_inc / I_coh = ", intensity_ratio

  ! Correlation functino given by Tr[ sigmap x exp(L*tau) x (sigmap x rho_ss) ]
  ! simgmam x rho_ss
  sigrho = MATMUL(rho_ss, sigmap)

  PRINT*, "Computing first order correlation"

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
    WRITE(1,*) k*dt, REAL(trace), IMAG(trace)

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

END PROGRAM first_order_correlation
