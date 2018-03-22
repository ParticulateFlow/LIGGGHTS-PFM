#ifndef LMP_TRI_MESH_CONTACTS_H
#define LMP_TRI_MESH_CONTACTS_H

#include <set>

namespace LAMMPS_NS
{
  class TriMesh;

  class TriMeshContacts
  {
   public:
    TriMeshContacts() : mesh_(NULL) {}
    TriMeshContacts(TriMesh * mesh) : mesh_(mesh) {}

    TriMesh * mesh_;
    std::set<int> contacts;
  };
}
#endif
