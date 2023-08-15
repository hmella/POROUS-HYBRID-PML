# Copyright (C) 2016-2022 by the multiphenics authors
#
# This file is part of multiphenics.
#
# multiphenics is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# multiphenics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with multiphenics. If not, see <http://www.gnu.org/licenses/>.
#

from dolfin import *
from multiphenics import *

# Helper function to generate subdomain restriction based on a gmsh subdomain id
def generate_subdomain_restriction(mesh, subdomains, subdomain_id):
    D = mesh.topology().dim()
    # Initialize empty restriction
    restriction = MeshRestriction(mesh, None)
    for d in range(D + 1):
        mesh_function_d = MeshFunction("bool", mesh, d)
        mesh_function_d.set_all(False)
        restriction.append(mesh_function_d)
    # Mark restriction mesh functions based on subdomain id
    for c in cells(mesh):
        if subdomains[c] == subdomain_id:
            restriction[D][c] = True
            for d in range(D):
                for e in entities(c, d):
                    restriction[d][e] = True
    # Return
    return restriction

# Helper function to generate interface restriction based on a pair of gmsh subdomain ids
def generate_interface_restriction(mesh, subdomains, subdomain_ids):
    assert isinstance(subdomain_ids, set)
    assert len(subdomain_ids) == 2
    D = mesh.topology().dim()
    # Initialize empty restriction
    restriction = MeshRestriction(mesh, None)
    for d in range(D + 1):
        mesh_function_d = MeshFunction("bool", mesh, d)
        mesh_function_d.set_all(False)
        restriction.append(mesh_function_d)
    # Mark restriction mesh functions based on subdomain ids (except the mesh function corresponding to dimension D, as it is trivially false)
    for f in facets(mesh):
        subdomains_ids_f = set(subdomains[c] for c in cells(f))
        assert len(subdomains_ids_f) in (1, 2)
        if subdomains_ids_f == subdomain_ids:
            restriction[D - 1][f] = True
            for d in range(D - 1):
                for e in entities(f, d):
                    restriction[d][e] = True
    # Return
    return restriction

if __name__=='__main__':

    # MPI parameters
    comm = MPI.comm_world
    rank = MPI.rank(comm)

    # Import mesh
    mesh = Mesh(comm)
    with XDMFFile(comm,"mesh/mesh_hybrid.xdmf") as infile:
        infile.read(mesh)

    # Markers for subdomains
    mvc = MeshValueCollection("size_t", mesh, mesh.geometry().dim())
    with XDMFFile(comm,"mesh/mf_hybrid.xdmf") as infile:
        infile.read(mvc, "markers")
    subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    # Generate restriction corresponding to interior subdomain (id = 1)
    RD_restriction = generate_subdomain_restriction(mesh, subdomains, 1)

    # Generate restriction corresponding to interior subdomain (id = 2)
    PML_restriction = generate_subdomain_restriction(mesh, subdomains, 2)

    # Generate restriction corresponding to interface between the two subdomains
    interface_restriction = generate_interface_restriction(mesh, subdomains, {1, 2})

    # Write out for simulation import (.xml) and visualization (.xdmf)
    File("mesh/mesh_RD_restriction.rtc.xml") << RD_restriction
    File("mesh/mesh_PML_restriction.rtc.xml") << PML_restriction
    File("mesh/mesh_interface_restriction.rtc.xml") << interface_restriction
    XDMFFile("mesh/mesh_RD_restriction.rtc.xdmf").write(RD_restriction)
    XDMFFile("mesh/mesh_PML_restriction.rtc.xdmf").write(PML_restriction)
    XDMFFile("mesh/mesh_interface_restriction.rtc.xdmf").write(interface_restriction)