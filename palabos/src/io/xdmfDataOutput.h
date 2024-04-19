/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef XDMF_DATA_OUTPUT_H
#define XDMF_DATA_OUTPUT_H

#ifdef HDF5

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "core/serializer.h"
#include "io/hdfWrapper.h"
#include "io/multiBlockWriter3D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

class ParallelXdmfDataWriter3D {
public:
    ParallelXdmfDataWriter3D(const std::string &fname);
    ~ParallelXdmfDataWriter3D();

    template <typename T>
    void writeDataField(MultiBlock3D &multiBlock, const std::string &field_name);

private:
    std::string xdmf_fname;
    std::string h5_fname;
    std::ofstream *fhandle;
    int field = 0;
};

}  // namespace plb

#endif

#endif
