#ifndef FILE_VTKOUTPUTQUAD_HPP
#define FILE_VTKOUTPUTQUAD_HPP

#include <comp.hpp>

namespace ngcomp
{
  template <int D>
  class VTKOutputQuad : public BaseVTKOutput
  {
  protected:
    enum CellType
    {
      VTK_VERTEX           = 1,
      VTK_LINE             = 3,
      VTK_TRIANGLE         = 5,
      VTK_QUAD             = 9,
      VTK_TETRA            = 10,
      VTK_HEXAHEDRON       = 12,
      VTK_WEDGE            = 13,
      VTK_PYRAMID          = 14,
      VTK_PENTAGONAL_PRISM = 15,
      VTK_HEXAGONAL_PRISM  = 16
    };

    shared_ptr<MeshAccess> ma = nullptr;
    Array<shared_ptr<CoefficientFunction>> coefs;
    Array<string> fieldnames;
    string filename;
    string grid_str;
    int subdivision;
    int only_element = -1;

    Array<shared_ptr<ValueField>> value_field;

    int output_cnt = 0;

    shared_ptr<ofstream> fileout;

  public:

    VTKOutputQuad(const Array<shared_ptr<CoefficientFunction>> &,
               const Flags &,shared_ptr<MeshAccess>);

    VTKOutputQuad(shared_ptr<MeshAccess>, const Array<shared_ptr<CoefficientFunction>> &,
               const Array<string> &, string, int, int);

    static int ElementTypeToVTKType(int et);
    void BuildGridString();
    void FillReferenceData(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_trigs);    
    // void FillReferenceData3D(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_tets);
    void PrintFieldData();

    virtual void Do (LocalHeap & lh, const BitArray * drawelems = 0);
  };
}

#endif
