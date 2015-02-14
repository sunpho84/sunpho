#ifndef _GEOMETRY_EO_H
#define _GEOMETRY_EO_H

#include "new_types/new_types_definitions.hpp"

namespace bissa
{
  int glblx_parity(int glx);
  int glb_coord_parity(coords c);
  void initialize_eo_bord_receivers_of_kind(MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base);
  void initialize_eo_bord_senders_of_kind(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *base);
  void set_eo_bord_senders_and_receivers(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base);
  void set_eo_geometry();
  void unset_eo_geometry();
}

#endif
