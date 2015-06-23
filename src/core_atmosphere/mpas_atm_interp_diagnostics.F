! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atm_interp_diagnostics
 use mpas_dmpar
 use mpas_kind_types
 use mpas_grid_types
 use mpas_constants

 implicit none
 private
 public:: interp_diagnostics

 contains

!==================================================================================================
 subroutine interp_diagnostics(mesh, state, time_lev, diag, diag_physics)
!==================================================================================================

!input arguments:
 type (mpas_pool_type), intent(in)  :: mesh
 type (mpas_pool_type), intent(in) :: state
 integer, intent(in) :: time_lev              ! which time level to use from state

!inout arguments:
 type (mpas_pool_type), intent(inout) :: diag
 type (mpas_pool_type), intent(inout) :: diag_physics

!local variables:
 integer :: iCell,iVert,iVertD,k,kk
 integer, pointer :: nCells, nVertLevels, nVertices, vertexDegree
 integer :: nVertLevelsP1
 integer, pointer :: index_qv
 integer, dimension(:,:), pointer :: cellsOnVertex

 type (field2DReal), pointer:: pressure_p_field

 real (kind=RKIND), dimension(:), pointer :: areaTriangle
 real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
 
 real (kind=RKIND), dimension(:,:), pointer :: exner, height
 real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p 
 real (kind=RKIND), dimension(:,:), pointer :: qvapor, relhum, theta_m, vorticity
 real (kind=RKIND), dimension(:,:), pointer :: umeridional, uzonal, vvel
 real (kind=RKIND), dimension(:,:,:), pointer :: scalars

 real (kind=RKIND), dimension(:), pointer :: temperature_020hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_025hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_050hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_100hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_150hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_200hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_300hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_400hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_500hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_600hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_700hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_850hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_925hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_1000hPa

 real (kind=RKIND), dimension(:), pointer :: qvapor_020hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_025hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_050hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_100hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_150hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_200hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_300hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_400hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_500hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_600hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_700hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_850hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_925hPa
 real (kind=RKIND), dimension(:), pointer :: qvapor_1000hPa

 real (kind=RKIND), dimension(:), pointer :: relhum_020hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_025hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_050hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_100hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_150hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_200hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_300hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_400hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_500hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_600hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_700hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_850hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_925hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_1000hPa

 real (kind=RKIND), dimension(:), pointer :: uzonal_020hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_025hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_050hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_100hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_150hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_200hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_300hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_400hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_500hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_600hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_700hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_850hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_925hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_1000hPa

 real (kind=RKIND), dimension(:), pointer :: umeridional_020hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_025hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_050hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_100hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_150hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_200hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_300hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_400hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_500hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_600hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_700hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_850hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_925hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_1000hPa

 real (kind=RKIND), dimension(:), pointer :: height_020hPa
 real (kind=RKIND), dimension(:), pointer :: height_025hPa
 real (kind=RKIND), dimension(:), pointer :: height_050hPa
 real (kind=RKIND), dimension(:), pointer :: height_100hPa
 real (kind=RKIND), dimension(:), pointer :: height_150hPa
 real (kind=RKIND), dimension(:), pointer :: height_200hPa
 real (kind=RKIND), dimension(:), pointer :: height_300hPa
 real (kind=RKIND), dimension(:), pointer :: height_400hPa
 real (kind=RKIND), dimension(:), pointer :: height_500hPa
 real (kind=RKIND), dimension(:), pointer :: height_600hPa
 real (kind=RKIND), dimension(:), pointer :: height_700hPa
 real (kind=RKIND), dimension(:), pointer :: height_850hPa
 real (kind=RKIND), dimension(:), pointer :: height_925hPa
 real (kind=RKIND), dimension(:), pointer :: height_1000hPa

 real (kind=RKIND), dimension(:), pointer :: w_020hPa
 real (kind=RKIND), dimension(:), pointer :: w_025hPa
 real (kind=RKIND), dimension(:), pointer :: w_050hPa
 real (kind=RKIND), dimension(:), pointer :: w_100hPa
 real (kind=RKIND), dimension(:), pointer :: w_150hPa
 real (kind=RKIND), dimension(:), pointer :: w_200hPa
 real (kind=RKIND), dimension(:), pointer :: w_300hPa
 real (kind=RKIND), dimension(:), pointer :: w_400hPa
 real (kind=RKIND), dimension(:), pointer :: w_500hPa
 real (kind=RKIND), dimension(:), pointer :: w_600hPa
 real (kind=RKIND), dimension(:), pointer :: w_700hPa
 real (kind=RKIND), dimension(:), pointer :: w_850hPa
 real (kind=RKIND), dimension(:), pointer :: w_925hPa
 real (kind=RKIND), dimension(:), pointer :: w_1000hPa

 real (kind=RKIND), dimension(:), pointer :: vorticity_200hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_500hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_700hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_850hPa

 real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressureCp1, pressure2, pressure_v, temperature

!local interpolated fields:
 integer:: nIntP
 real(kind=RKIND):: w1,w2,z0,z1,z2
 real(kind=RKIND),dimension(:,:),allocatable:: field_in,press_in, temp_in, tempp1_in
 real(kind=RKIND),dimension(:,:),allocatable:: field_interp,press_interp
 real(kind=RKIND),dimension(:),allocatable:: tmp
 real(kind=RKIND):: dum

!--------------------------------------------------------------------------------------------------

! write(0,*)
! write(0,*) '--- enter subroutine interp_diagnostics:'

 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
 call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
 call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
 call mpas_pool_get_dimension(state, 'index_qv', index_qv)

 nVertLevelsP1 = nVertLevels + 1

 call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
 call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
 call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)

 call mpas_pool_get_array(mesh, 'zgrid', height)
 call mpas_pool_get_array(state, 'w', vvel, time_lev)
 call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
 call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

 qvapor  => scalars(index_qv,:,:)    !MGD does this actually work?
 
 call mpas_pool_get_field(diag, 'pressure_p', pressure_p_field)
 call mpas_dmpar_exch_halo_field(pressure_p_field)

 call mpas_pool_get_array(diag, 'exner', exner)
 call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
 call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
 call mpas_pool_get_array(diag, 'vorticity', vorticity)
 call mpas_pool_get_array(diag, 'uReconstructMeridional', umeridional)
 call mpas_pool_get_array(diag, 'uReconstructZonal', uzonal)
 call mpas_pool_get_array(diag_physics, 'relhum', relhum)

 call mpas_pool_get_array(diag, 'temperature_020hPa', temperature_020hPa)
 call mpas_pool_get_array(diag, 'temperature_025hPa', temperature_025hPa)
 call mpas_pool_get_array(diag, 'temperature_050hPa', temperature_050hPa)
 call mpas_pool_get_array(diag, 'temperature_100hPa', temperature_100hPa)
 call mpas_pool_get_array(diag, 'temperature_200hPa', temperature_200hPa)
 call mpas_pool_get_array(diag, 'temperature_300hPa', temperature_300hPa)
 call mpas_pool_get_array(diag, 'temperature_400hPa', temperature_400hPa)
 call mpas_pool_get_array(diag, 'temperature_500hPa', temperature_500hPa)
 call mpas_pool_get_array(diag, 'temperature_600hPa', temperature_600hPa)
 call mpas_pool_get_array(diag, 'temperature_700hPa', temperature_700hPa)
 call mpas_pool_get_array(diag, 'temperature_850hPa', temperature_850hPa)
 call mpas_pool_get_array(diag, 'temperature_925hPa', temperature_925hPa)
 call mpas_pool_get_array(diag, 'temperature_1000hPa', temperature_1000hPa)

 call mpas_pool_get_array(diag, 'relhum_020hPa', relhum_020hPa)
 call mpas_pool_get_array(diag, 'relhum_025hPa', relhum_025hPa)
 call mpas_pool_get_array(diag, 'relhum_050hPa', relhum_050hPa)
 call mpas_pool_get_array(diag, 'relhum_100hPa', relhum_100hPa)
 call mpas_pool_get_array(diag, 'relhum_200hPa', relhum_200hPa)
 call mpas_pool_get_array(diag, 'relhum_300hPa', relhum_300hPa)
 call mpas_pool_get_array(diag, 'relhum_400hPa', relhum_400hPa)
 call mpas_pool_get_array(diag, 'relhum_500hPa', relhum_500hPa)
 call mpas_pool_get_array(diag, 'relhum_600hPa', relhum_600hPa)
 call mpas_pool_get_array(diag, 'relhum_700hPa', relhum_700hPa)
 call mpas_pool_get_array(diag, 'relhum_850hPa', relhum_850hPa)
 call mpas_pool_get_array(diag, 'relhum_925hPa', relhum_925hPa)
 call mpas_pool_get_array(diag, 'relhum_1000hPa', relhum_1000hPa)

 call mpas_pool_get_array(diag, 'qvapor_020hPa', qvapor_020hPa)
 call mpas_pool_get_array(diag, 'qvapor_025hPa', qvapor_025hPa)
 call mpas_pool_get_array(diag, 'qvapor_050hPa', qvapor_050hPa)
 call mpas_pool_get_array(diag, 'qvapor_100hPa', qvapor_100hPa)
 call mpas_pool_get_array(diag, 'qvapor_200hPa', qvapor_200hPa)
 call mpas_pool_get_array(diag, 'qvapor_300hPa', qvapor_300hPa)
 call mpas_pool_get_array(diag, 'qvapor_400hPa', qvapor_400hPa)
 call mpas_pool_get_array(diag, 'qvapor_500hPa', qvapor_500hPa)
 call mpas_pool_get_array(diag, 'qvapor_600hPa', qvapor_600hPa)
 call mpas_pool_get_array(diag, 'qvapor_700hPa', qvapor_700hPa)
 call mpas_pool_get_array(diag, 'qvapor_850hPa', qvapor_850hPa)
 call mpas_pool_get_array(diag, 'qvapor_925hPa', qvapor_925hPa)
 call mpas_pool_get_array(diag, 'qvapor_1000hPa', qvapor_1000hPa)

 call mpas_pool_get_array(diag, 'uzonal_020hPa', uzonal_020hPa)
 call mpas_pool_get_array(diag, 'uzonal_025hPa', uzonal_025hPa)
 call mpas_pool_get_array(diag, 'uzonal_050hPa', uzonal_050hPa)
 call mpas_pool_get_array(diag, 'uzonal_100hPa', uzonal_100hPa)
 call mpas_pool_get_array(diag, 'uzonal_200hPa', uzonal_200hPa)
 call mpas_pool_get_array(diag, 'uzonal_300hPa', uzonal_300hPa)
 call mpas_pool_get_array(diag, 'uzonal_400hPa', uzonal_400hPa)
 call mpas_pool_get_array(diag, 'uzonal_500hPa', uzonal_500hPa)
 call mpas_pool_get_array(diag, 'uzonal_600hPa', uzonal_600hPa)
 call mpas_pool_get_array(diag, 'uzonal_700hPa', uzonal_700hPa)
 call mpas_pool_get_array(diag, 'uzonal_850hPa', uzonal_850hPa)
 call mpas_pool_get_array(diag, 'uzonal_925hPa', uzonal_925hPa)
 call mpas_pool_get_array(diag, 'uzonal_1000hPa', uzonal_1000hPa)

 call mpas_pool_get_array(diag, 'umeridional_020hPa', umeridional_020hPa)
 call mpas_pool_get_array(diag, 'umeridional_025hPa', umeridional_025hPa)
 call mpas_pool_get_array(diag, 'umeridional_050hPa', umeridional_050hPa)
 call mpas_pool_get_array(diag, 'umeridional_100hPa', umeridional_100hPa)
 call mpas_pool_get_array(diag, 'umeridional_200hPa', umeridional_200hPa)
 call mpas_pool_get_array(diag, 'umeridional_300hPa', umeridional_300hPa)
 call mpas_pool_get_array(diag, 'umeridional_400hPa', umeridional_400hPa)
 call mpas_pool_get_array(diag, 'umeridional_500hPa', umeridional_500hPa)
 call mpas_pool_get_array(diag, 'umeridional_600hPa', umeridional_600hPa)
 call mpas_pool_get_array(diag, 'umeridional_700hPa', umeridional_700hPa)
 call mpas_pool_get_array(diag, 'umeridional_850hPa', umeridional_850hPa)
 call mpas_pool_get_array(diag, 'umeridional_925hPa', umeridional_925hPa)
 call mpas_pool_get_array(diag, 'umeridional_1000hPa', umeridional_1000hPa)

 call mpas_pool_get_array(diag, 'height_020hPa', height_020hPa)
 call mpas_pool_get_array(diag, 'height_025hPa', height_025hPa)
 call mpas_pool_get_array(diag, 'height_050hPa', height_050hPa)
 call mpas_pool_get_array(diag, 'height_100hPa', height_100hPa)
 call mpas_pool_get_array(diag, 'height_200hPa', height_200hPa)
 call mpas_pool_get_array(diag, 'height_300hPa', height_300hPa)
 call mpas_pool_get_array(diag, 'height_400hPa', height_400hPa)
 call mpas_pool_get_array(diag, 'height_500hPa', height_500hPa)
 call mpas_pool_get_array(diag, 'height_600hPa', height_600hPa)
 call mpas_pool_get_array(diag, 'height_700hPa', height_700hPa)
 call mpas_pool_get_array(diag, 'height_850hPa', height_850hPa)
 call mpas_pool_get_array(diag, 'height_925hPa', height_925hPa)
 call mpas_pool_get_array(diag, 'height_1000hPa', height_1000hPa)

 call mpas_pool_get_array(diag, 'w_020hPa', w_020hPa)
 call mpas_pool_get_array(diag, 'w_025hPa', w_025hPa)
 call mpas_pool_get_array(diag, 'w_050hPa', w_050hPa)
 call mpas_pool_get_array(diag, 'w_100hPa', w_100hPa)
 call mpas_pool_get_array(diag, 'w_200hPa', w_200hPa)
 call mpas_pool_get_array(diag, 'w_300hPa', w_300hPa)
 call mpas_pool_get_array(diag, 'w_400hPa', w_400hPa)
 call mpas_pool_get_array(diag, 'w_500hPa', w_500hPa)
 call mpas_pool_get_array(diag, 'w_600hPa', w_600hPa)
 call mpas_pool_get_array(diag, 'w_700hPa', w_700hPa)
 call mpas_pool_get_array(diag, 'w_850hPa', w_850hPa)
 call mpas_pool_get_array(diag, 'w_925hPa', w_925hPa)
 call mpas_pool_get_array(diag, 'w_1000hPa', w_1000hPa)

 call mpas_pool_get_array(diag, 'vorticity_200hPa', vorticity_200hPa)
 call mpas_pool_get_array(diag, 'vorticity_500hPa', vorticity_500hPa)
 call mpas_pool_get_array(diag, 'vorticity_700hPa', vorticity_700hPa)
 call mpas_pool_get_array(diag, 'vorticity_850hPa', vorticity_850hPa)

 if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCells)      )
 if(.not.allocated(pressureCp1) ) allocate(pressureCp1(nVertLevels,nCells+1) )
 if(.not.allocated(pressure2)   ) allocate(pressure2(nVertLevelsP1,nCells)   )
 if(.not.allocated(pressure_v)  ) allocate(pressure_v(nVertLevels,nVertices) )
 if(.not.allocated(temperature) ) allocate(temperature(nVertLevels,nCells)   )

!calculation of total pressure at cell centers (at mass points):
 do iCell = 1, nCells
 do k = 1, nVertLevels
    pressure(k,iCell)    = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
    pressureCp1(k,iCell) = pressure(k,iCell)
 enddo
 enddo
 do iCell = nCells+1, nCells+1
 do k = 1, nVertLevels
    pressureCp1(k,iCell)   = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
 enddo
 enddo

!calculation of total pressure at cell centers (at vertical velocity points):
 k = nVertLevelsP1
 do iCell = 1, nCells
    z0 = height(k,iCell)
    z1 = 0.5*(height(k,iCell)+height(k-1,iCell)) 
    z2 = 0.5*(height(k-1,iCell)+height(k-2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pressure2(k,iCell) = exp(w1*log(pressure(k-1,iCell))+w2*log(pressure(k-2,iCell)))
 enddo
 do k = 2, nVertLevels
 do iCell = 1, nCells
    w1 = (height(k,iCell)-height(k-1,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
    w2 = (height(k+1,iCell)-height(k,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
    pressure2(k,iCell) = w1*pressure(k,iCell) + w2*pressure(k-1,iCell)
 enddo
 enddo
 k = 1
 do iCell = 1, nCells
    z0 = height(k,iCell)
    z1 = 0.5*(height(k,iCell)+height(k+1,iCell)) 
    z2 = 0.5*(height(k+1,iCell)+height(k+2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    pressure2(k,iCell) = w1*pressure(k,iCell)+w2*pressure(k+1,iCell)
 enddo

!calculation of total pressure at cell vertices (at mass points):
 do iVert = 1, nVertices
    pressure_v(:,iVert) = 0._RKIND

    do k = 1, nVertLevels
    do iVertD = 1, vertexDegree
       pressure_v(k,iVert) = pressure_v(k,iVert) &
               + kiteAreasOnVertex(iVertD,iVert)*pressureCp1(k,cellsOnVertex(iVertD,iVert))
    enddo
    pressure_v(k,iVert) = pressure_v(k,iVert) / areaTriangle(iVert)
    enddo
 enddo

!calculation of temperature at cell centers:
 do iCell = 1,nCells
 do k = 1,nVertLevels
    temperature(k,iCell) = (theta_m(k,iCell)/(1._RKIND+rvord*qvapor(k,iCell)))*exner(k,iCell) 
 enddo
 enddo

!interpolation to fixed pressure levels for fields located at cells centers and at mass points:
 nIntP = 14
 if(.not.allocated(field_interp)) allocate(field_interp(nCells,nIntP) )
 if(.not.allocated(press_interp)) allocate(press_interp(nCells,nIntP) )
 do iCell = 1, nCells
    press_interp(iCell,1) = 020.0_RKIND
    press_interp(iCell,2) = 025.0_RKIND
    press_interp(iCell,3) = 050.0_RKIND
    press_interp(iCell,4) = 100.0_RKIND
    press_interp(iCell,5) = 150.0_RKIND
    press_interp(iCell,6) = 200.0_RKIND
    press_interp(iCell,7) = 300.0_RKIND
    press_interp(iCell,8) = 400.0_RKIND
    press_interp(iCell,9) = 500.0_RKIND
    press_interp(iCell,10) = 600.0_RKIND
    press_interp(iCell,11) = 700.0_RKIND
    press_interp(iCell,12) = 850.0_RKIND
    press_interp(iCell,13) = 925.0_RKIND
    press_interp(iCell,14) = 1000.0_RKIND
 enddo

 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    press_in(iCell,kk) = pressure(k,iCell)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
!... temperature:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = temperature(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure_fim(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp,temp_in,1,0) 
 temperature_020hPa(1:nCells) = field_interp(1:nCells,1)
 temperature_025hPa(1:nCells) = field_interp(1:nCells,2)
 temperature_050hPa(1:nCells) = field_interp(1:nCells,3)
 temperature_100hPa(1:nCells) = field_interp(1:nCells,4)
 temperature_150hPa(1:nCells) = field_interp(1:nCells,5)
 temperature_200hPa(1:nCells) = field_interp(1:nCells,6)
 temperature_300hPa(1:nCells) = field_interp(1:nCells,7)
 temperature_400hPa(1:nCells) = field_interp(1:nCells,8)
 temperature_500hPa(1:nCells) = field_interp(1:nCells,9)
 temperature_600hPa(1:nCells) = field_interp(1:nCells,10)
 temperature_700hPa(1:nCells) = field_interp(1:nCells,11)
 temperature_850hPa(1:nCells) = field_interp(1:nCells,12)
 temperature_925hPa(1:nCells) = field_interp(1:nCells,13)
 temperature_1000hPa(1:nCells)= field_interp(1:nCells,14)
! write(0,*) '--- end interpolate temperature:'

!... mixing Ratio
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = qvapor(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 qvapor_020hPa(1:nCells) = field_interp(1:nCells,1)
 qvapor_025hPa(1:nCells) = field_interp(1:nCells,2)
 qvapor_050hPa(1:nCells) = field_interp(1:nCells,3)
 qvapor_100hPa(1:nCells) = field_interp(1:nCells,4)
 qvapor_150hPa(1:nCells) = field_interp(1:nCells,5)
 qvapor_200hPa(1:nCells) = field_interp(1:nCells,6)
 qvapor_300hPa(1:nCells) = field_interp(1:nCells,7)
 qvapor_400hPa(1:nCells) = field_interp(1:nCells,8)
 qvapor_500hPa(1:nCells) = field_interp(1:nCells,9)
 qvapor_600hPa(1:nCells) = field_interp(1:nCells,10)
 qvapor_700hPa(1:nCells) = field_interp(1:nCells,11)
 qvapor_850hPa(1:nCells) = field_interp(1:nCells,12)
 qvapor_925hPa(1:nCells) = field_interp(1:nCells,13)
 qvapor_1000hPa(1:nCells)= field_interp(1:nCells,14)
 write(0,*) '--- end interpolate mixing ratio:'

!... relative humidity:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = relhum(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 relhum_020hPa(1:nCells) = field_interp(1:nCells,1)
 relhum_025hPa(1:nCells) = field_interp(1:nCells,2)
 relhum_050hPa(1:nCells) = field_interp(1:nCells,3)
 relhum_100hPa(1:nCells) = field_interp(1:nCells,4)
 relhum_150hPa(1:nCells) = field_interp(1:nCells,5)
 relhum_200hPa(1:nCells) = field_interp(1:nCells,6)
 relhum_300hPa(1:nCells) = field_interp(1:nCells,7)
 relhum_400hPa(1:nCells) = field_interp(1:nCells,8)
 relhum_500hPa(1:nCells) = field_interp(1:nCells,9)
 relhum_600hPa(1:nCells) = field_interp(1:nCells,10)
 relhum_700hPa(1:nCells) = field_interp(1:nCells,11)
 relhum_850hPa(1:nCells) = field_interp(1:nCells,12)
 relhum_925hPa(1:nCells) = field_interp(1:nCells,13)
 relhum_1000hPa(1:nCells)= field_interp(1:nCells,14)
! write(0,*) '--- end interpolate relative humidity:'

!... u zonal wind:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = uzonal(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 uzonal_020hPa(1:nCells) = field_interp(1:nCells,1)
 uzonal_025hPa(1:nCells) = field_interp(1:nCells,2)
 uzonal_050hPa(1:nCells) = field_interp(1:nCells,3)
 uzonal_100hPa(1:nCells) = field_interp(1:nCells,4)
 uzonal_150hPa(1:nCells) = field_interp(1:nCells,5)
 uzonal_200hPa(1:nCells) = field_interp(1:nCells,6)
 uzonal_300hPa(1:nCells) = field_interp(1:nCells,7)
 uzonal_400hPa(1:nCells) = field_interp(1:nCells,8)
 uzonal_500hPa(1:nCells) = field_interp(1:nCells,9)
 uzonal_600hPa(1:nCells) = field_interp(1:nCells,10)
 uzonal_700hPa(1:nCells) = field_interp(1:nCells,11)
 uzonal_850hPa(1:nCells) = field_interp(1:nCells,12)
 uzonal_925hPa(1:nCells) = field_interp(1:nCells,13)
 uzonal_1000hPa(1:nCells)= field_interp(1:nCells,14)
! write(0,*) '--- end interpolate zonal wind:'

!... u meridional wind:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = umeridional(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 umeridional_020hPa(1:nCells) = field_interp(1:nCells,1)
 umeridional_025hPa(1:nCells) = field_interp(1:nCells,2)
 umeridional_050hPa(1:nCells) = field_interp(1:nCells,3)
 umeridional_100hPa(1:nCells) = field_interp(1:nCells,4)
 umeridional_150hPa(1:nCells) = field_interp(1:nCells,5)
 umeridional_200hPa(1:nCells) = field_interp(1:nCells,6)
 umeridional_300hPa(1:nCells) = field_interp(1:nCells,7)
 umeridional_400hPa(1:nCells) = field_interp(1:nCells,8)
 umeridional_500hPa(1:nCells) = field_interp(1:nCells,9)
 umeridional_600hPa(1:nCells) = field_interp(1:nCells,10)
 umeridional_700hPa(1:nCells) = field_interp(1:nCells,11)
 umeridional_850hPa(1:nCells) = field_interp(1:nCells,12)
 umeridional_925hPa(1:nCells) = field_interp(1:nCells,13)
 umeridional_1000hPa(1:nCells)= field_interp(1:nCells,14)
! write(0,*) '--- end interpolate meridional wind:'

 if(allocated(field_in)) deallocate(field_in)
 if(allocated(press_in)) deallocate(press_in)

!interpolation to fixed pressure levels for fields located at cells centers and at vertical
!velocity points:
 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevelsP1))
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    press_in(iCell,kk) = pressure2(k,iCell)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevelsP1))
 !... height:
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    field_in(iCell,kk) = height(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
 call interp_tofixed_pressure_fim(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp,temp_in,0,1) 
 height_020hPa(1:nCells) = field_interp(1:nCells,1)
 height_025hPa(1:nCells) = field_interp(1:nCells,2)
 height_050hPa(1:nCells) = field_interp(1:nCells,3)
 height_100hPa(1:nCells) = field_interp(1:nCells,4)
 height_150hPa(1:nCells) = field_interp(1:nCells,5)
 height_200hPa(1:nCells) = field_interp(1:nCells,6)
 height_300hPa(1:nCells) = field_interp(1:nCells,7)
 height_400hPa(1:nCells) = field_interp(1:nCells,8)
 height_500hPa(1:nCells) = field_interp(1:nCells,9)
 height_600hPa(1:nCells) = field_interp(1:nCells,10)
 height_700hPa(1:nCells) = field_interp(1:nCells,11)
 height_850hPa(1:nCells) = field_interp(1:nCells,12)
 height_925hPa(1:nCells) = field_interp(1:nCells,13)
 height_1000hPa(1:nCells)= field_interp(1:nCells,14)
! write(0,*) '--- end interpolate height:'
 
!... vertical velocity
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    field_in(iCell,kk) = vvel(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
 w_020hPa(1:nCells) = field_interp(1:nCells,1)
 w_025hPa(1:nCells) = field_interp(1:nCells,2)
 w_050hPa(1:nCells) = field_interp(1:nCells,3)
 w_100hPa(1:nCells) = field_interp(1:nCells,4)
 w_150hPa(1:nCells) = field_interp(1:nCells,5)
 w_200hPa(1:nCells) = field_interp(1:nCells,6)
 w_300hPa(1:nCells) = field_interp(1:nCells,7)
 w_400hPa(1:nCells) = field_interp(1:nCells,8)
 w_500hPa(1:nCells) = field_interp(1:nCells,9)
 w_600hPa(1:nCells) = field_interp(1:nCells,10)
 w_700hPa(1:nCells) = field_interp(1:nCells,11)
 w_850hPa(1:nCells) = field_interp(1:nCells,12)
 w_925hPa(1:nCells) = field_interp(1:nCells,13)
 w_1000hPa(1:nCells)= field_interp(1:nCells,14)
! write(0,*) '--- end interpolate vertical velocity:'
 
 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)

!interpolation to fixed pressure levels for fields located at cell vertices and at mass points:
 nIntP = 4
 if(.not.allocated(field_interp)) allocate(field_interp(nVertices,nIntP) )
 if(.not.allocated(press_interp)) allocate(press_interp(nVertices,nIntP) )
 do iVert = 1, nVertices
    press_interp(iVert,1) = 200.0_RKIND
    press_interp(iVert,2) = 500.0_RKIND
    press_interp(iVert,3) = 700.0_RKIND
    press_interp(iVert,4) = 850.0_RKIND
 enddo

 if(allocated(field_in)) deallocate(field_in)
 if(allocated(press_in)) deallocate(press_in)

 if(.not.allocated(press_in)) allocate(press_in(nVertices,nVertLevels))
 do iVert = 1, nVertices
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    press_in(iVert,kk) = pressure_v(k,iVert)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nVertices,nVertLevels))
!... relative vorticity:
 do iVert = 1, nVertices
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iVert,kk) = vorticity(k,iVert)
 enddo
 enddo
 call interp_tofixed_pressure(nVertices,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 vorticity_200hPa(1:nVertices) = field_interp(1:nVertices,1)
 vorticity_500hPa(1:nVertices) = field_interp(1:nVertices,2)
 vorticity_700hPa(1:nVertices) = field_interp(1:nVertices,3)
 vorticity_850hPa(1:nVertices) = field_interp(1:nVertices,4)
! write(0,*) '--- end interpolate relative vorticity:'

 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)
 if(allocated(pressure)    ) deallocate(pressure    )
 if(allocated(pressureCp1) ) deallocate(pressureCp1 )
 if(allocated(pressure2)   ) deallocate(pressure2   )
 if(allocated(pressure_v)  ) deallocate(pressure_v  )
 if(allocated(temperature) ) deallocate(temperature )

!formats:
! 201 format(i5,4(1x,e15.8))

 end subroutine interp_diagnostics

!==================================================================================================
 subroutine interp_tofixed_pressure_fim(ncol,nlev_in,nlev_out,p_in,field_in,p_out,field_out, &
             tinp,it,iz)
!==================================================================================================

   use mpas_atmphys_constants, only:rcp, p0,rdg

!input arguments:
 integer,intent(in):: ncol,nlev_in,nlev_out
 integer, intent(in) :: it,iz

 real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: p_in,field_in,tinp
 real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: p_out

!output arguments:
 real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out

!local variables:
 integer:: i1,i2
 integer:: icol,k,kk
 integer:: var,ka,kaa,kb
 integer:: kkstart,kount
 integer,dimension(ncol):: kupper

 real(kind=RKIND):: dpl,dpu,thav
 real(kind=RKIND),dimension(ncol,nlev_in) :: pres_in
 real(kind=RKIND),dimension(ncol,nlev_in) :: thinp
 real(kind=RKIND),dimension(ncol,nlev_out):: pres_out

 real(kind=RKIND):: pk1,pk2

!--------------------------------------------------------------------------------------------------

!formats:
 201 format(i5,8(1x,e15.8))

!write(0,*)
!write(0,*) '--- enter subroutine interp_tofixed_pressure_fim:'
!write(0,*) '... ncol     = ',ncol
!write(0,*) '... nlev_in  = ',nlev_in
!write(0,*) '... nlev_out = ',nlev_out
 i1=1 ; i2=ncol/2


 do icol = 1, ncol
    kupper(icol) = 1
    pres_in(icol,:) = (p_in(icol,:)/(0.010_RKIND*p0) ) ** rcp
    pres_out(icol,:) = (p_out(icol,:)/(0.010_RKIND*p0) ) ** rcp
    thinp(icol,:) = tinp(icol,:) / pres_in(icol,:)
 enddo

! do k = 1, nlev_in
!    write(0,201)
!    k,p_in(i1,k),field_in(i1,k),tinp(i1,k),pres_in(i1,k),thinp(i1,k)
! enddo
! write(0,*)
!
! do k = 1, nlev_out
!    write(0,201) k,p_out(i1,k),pres_out(i1,k)
! enddo
! write(0,*)

 do k = 1, nlev_out

    kkstart = nlev_in
    do icol = 1, ncol
       kkstart = min0(kkstart,kupper(icol))
    enddo
    kount = 0

    do kk = kkstart, nlev_in-1
       do icol = 1, ncol
          if(pres_out(icol,k).gt.pres_in(icol,kk).and.pres_out(icol,k).le.pres_in(icol,kk+1)) then
             kupper(icol) = kk
             kount = kount + 1
!            write(0,201)
!            kupper(icol),pres_out(icol,k),pres_in(icol,kk),pres_in(icol,kk+1)
          endif
       enddo

!      if(kount.eq.ncol) then
!         do icol = 1, ncol
!            dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
!            dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
!            write(0,201) icol,pres_out(icol,k),pres_in(icol,kupper(icol)),pres_in(icol,kupper(icol)+1)
!            field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
!                              + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
!         end do
!         goto 35
!       end if
    enddo

    do icol = 1, ncol
       if(pres_out(icol,k) .lt. pres_in(icol,1)) then
          field_out(icol,k) = field_in(icol,1)*pres_out(icol,k)/pres_in(icol,1)

       elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in) .AND. (it+iz) .eq. 0 ) then
! makes the output field a constant below the surface!
          field_out(icol,k) = field_in(icol,nlev_in)
!      elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in) .AND. (it+iz) .gt. 0 )
!      then
       elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in) .AND. (it+iz) .gt. 0 ) then
! Extrapolate down below surface
!         write(0,*) ' Do interpolation
!         ',pres_out(icol,k),pres_in(icol,nlev_in),it,iz
!         ka=max(      1,kk-2)
!         kb=min(nlev_in,kk  )
          kaa=nlev_in-3
          ka=nlev_in-2
          kb=nlev_in
!         if(iz.eq.1) write(0,*) "ka  ",ka,kb,kk,nlev_in

! --- interpolate hydrostatically constrained pot.temperature
          var=.5*min(thinp(icol,kb)-thinp(icol,kb-1), thinp(icol,kb-1)-thinp(icol,ka))
!         write(0,*) "var ",var,thinp(icol,kk),thinp(icol,ka),thinp(icol,kb)

! --- assume linear variation within layer (as in PLM)
          pk1=pres_in(icol,ka  )
          pk2=pres_in(icol,kb  )
!          if(pk2 .lt. .95) then
!          write(0,*) "pk1 ",pk1,pk2,thinp(i1,kk)
!          write(0,*) "pk2",pk2,pres_out(icol,k),pk2-pres_out(icol,k)
!          write(0,*) "dp ",pres_out(icol,k),pk1,pres_out(icol,k)-pk1
!          endif

!         field_out(icol,k)=((thinp(icol,kk-1)-var)*(pres_out(icol,k)-pk2)
!         &
!                           +(thinp(icol,kk-1)+var)*(pk1-pres_out(icol,k)))/(pk1-pk2)
          field_out(icol,k)=((thinp(icol,nlev_in-1)-var)*(pres_out(icol,k)-pk2) &
                            +(thinp(icol,nlev_in-1)+var)*(pk1-pres_out(icol,k)))/(pk1-pk2)
!          if(pk2 .lt. .95) then
!         write(0,*) 'field out',field_out(icol,k)
!          endif

! --- interpolate geopotential
          if(iz .eq. 1) then
!           thav=thinp(icol,kk)                            ! traditional method
!           write(0,*) 'interpolated height'
!           thav = .5*(thinp(icol,kk-1)+var+field_out(icol,k))
            thav = .5*(thinp(icol,nlev_in-1)+var+field_out(icol,k))
            thav = rdg/rcp * thav
            field_out(icol,k)=field_in(icol,ka)-thav*(pres_out(icol,k) - pres_in(icol,ka))
!           write(0,'2(f15.5,f15.9)')
!           field_in(icol,kb),pres_in(icol,kb),field_out(icol,k),pres_out(icol,k)
          endif

! --- convert back to temperature from potential temperature
          if(it .eq. 1) then
            field_out(icol,k)=field_out(icol,k) * pres_out(icol,k)
            write(0,'(A20,6f9.3)') 'interp temp ',field_out(icol,k),tinp(icol,ka),thinp(icol,nlev_in-1), pres_out(icol,k),pk2,pk1
          endif

!         write(0,201) kb,field_in(i1,ka),field_in(i1,kb),field_out(i1,k)
!         write(0,*) ka,kb,var,pk1,pk2
!         write(0,*) thinp(i1,kk),thav,pres_in(i1,kk),pres_out(i1,ka)
!         write(0,*) '  '
       else
          if( pres_out(icol,k).LT.pres_in(icol,kupper(icol))  )  then
            write(*,*) ' ERROR: negative dpu ',pres_out(icol,k), pres_in(icol,kupper(icol)  )
            stop
          endif
          if( pres_out(icol,k).GT.pres_in(icol,kupper(icol)+1)) then
            write(*,*) ' ERROR: negative dpl ',pres_out(icol,k), pres_in(icol,kupper(icol)+1)
            stop
          endif

          dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
          dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
          field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                            + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
       endif
    enddo

 35 continue

!   write(0,201)
!   kupper(i1),pres_out(i1,k),pres_in(i1,kupper(i1)),pres_in(i1,kupper(i1)+1), &
!                field_out(i1,k),field_in(i1,kupper(i1)),field_in(i1,kupper(i1)+1)
!   write(0,201)
!   kupper(i2),pres_out(i2,k),pres_in(i2,kupper(i2)),pres_in(i2,kupper(i2)+1), &
!                field_out(i2,k),field_in(i2,kupper(i2)),field_in(i2,kupper(i2)+1)

 enddo

 end subroutine interp_tofixed_pressure_fim
!==================================================================================================
 subroutine interp_tofixed_pressure(ncol,nlev_in,nlev_out,pres_in,field_in,pres_out,field_out)
!==================================================================================================

!input arguments:
 integer,intent(in):: ncol,nlev_in,nlev_out

 real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: pres_in,field_in
 real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: pres_out

!output arguments:
 real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out

!local variables:
! integer:: i1,i2,icol,k,kk
 integer:: icol,k,kk
 integer:: kkstart,kount
 integer,dimension(ncol):: kupper

 real(kind=RKIND):: dpl,dpu

!--------------------------------------------------------------------------------------------------

!formats:
! 201 format(i5,8(1x,e15.8))

!write(0,*)
!write(0,*) '--- enter subroutine interp_tofixed_pressure:'
!write(0,*) '... ncol     = ',ncol
!write(0,*) '... nlev_in  = ',nlev_in
!write(0,*) '... nlev_out = ',nlev_out
!i1=1 ; i2=ncol
!do k = 1, nlev_in
!   write(0,201) k,pres_in(i1,k),field_in(i1,k),pres_in(i2,k),field_in(i2,k)
!enddo
!write(0,*)

 do icol = 1, ncol
    kupper(icol) = 1
 enddo

 do k = 1, nlev_out 

    kkstart = nlev_in
    do icol = 1, ncol
       kkstart = min0(kkstart,kupper(icol))
    enddo
    kount = 0

    do kk = kkstart, nlev_in-1
       do icol = 1, ncol
          if(pres_out(icol,k).gt.pres_in(icol,kk).and.pres_out(icol,k).le.pres_in(icol,kk+1)) then
             kupper(icol) = kk
             kount = kount + 1
!            write(0,201) kupper(icol),pres_out(icol,k),pres_in(icol,kk),pres_in(icol,kk+1)
          endif
       enddo

       if(kount.eq.ncol) then
          do icol = 1, ncol
             dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
             dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
             field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                               + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
          end do
          goto 35
        end if
    enddo

    do icol = 1, ncol
       if(pres_out(icol,k) .lt. pres_in(icol,1)) then
          field_out(icol,k) = field_in(icol,1)*pres_out(icol,k)/pres_in(icol,1)
       elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in)) then
          field_out(icol,k) = field_in(icol,nlev_in)
       else
          dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
          dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
          field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                            + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
       endif
    enddo

 35 continue
!   write(0,201) kupper(i1),pres_out(i1,k),pres_in(i1,kupper(i1)),pres_in(i1,kupper(i1)+1), &
!                field_out(i1,k),field_in(i1,kupper(i1)),field_in(i1,kupper(i1)+1)
!   write(0,201) kupper(i2),pres_out(i2,k),pres_in(i2,kupper(i2)),pres_in(i2,kupper(i2)+1), &
!                field_out(i2,k),field_in(i2,kupper(i2)),field_in(i2,kupper(i2)+1)

 enddo

 end subroutine interp_tofixed_pressure

!==================================================================================================
 end module mpas_atm_interp_diagnostics
!==================================================================================================
