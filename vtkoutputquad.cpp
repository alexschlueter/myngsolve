#include "vtkoutputquad.hpp"

namespace ngcomp
{
  template <int D>
  VTKOutputQuad<D>::VTKOutputQuad(const Array<shared_ptr<CoefficientFunction>> & a_coefs,
                           const Flags & flags,
                           shared_ptr<MeshAccess> ama )
    : VTKOutputQuad(ama, a_coefs,
                flags.GetStringListFlag ("fieldnames" ),
                flags.GetStringFlag ("filename","output"),
                (int) flags.GetNumFlag ( "subdivision", 0),
                (int) flags.GetNumFlag ( "only_element", -1))
  {;}

  template <int D>
  VTKOutputQuad<D>::VTKOutputQuad(shared_ptr<MeshAccess> ama,
                           const Array<shared_ptr<CoefficientFunction>> & a_coefs,
                           const Array<string> & a_field_names,
                           string a_filename, int a_subdivision, int a_only_element)
    : ma(ama), coefs(a_coefs), fieldnames(a_field_names),
      filename(a_filename), subdivision(a_subdivision), only_element(a_only_element)
  {
    BuildGridString();
    value_field.SetSize(a_coefs.Size());
    for (int i = 0; i < a_coefs.Size(); i++)
      if (fieldnames.Size() > i)
        value_field[i] = make_shared<ValueField>(coefs[i]->Dimension(),fieldnames[i]);
      else
        value_field[i] = make_shared<ValueField>(coefs[i]->Dimension(),"dummy" + to_string(i));
  }


  template <int D>
  void VTKOutputQuad<D>::BuildGridString()
  {
    ostringstream ss;
    // header:
    ss << "# vtk DataFile Version 3.0" << endl;
    ss << "vtk output" << endl;
    ss << "ASCII" << endl;
    ss << "DATASET UNSTRUCTURED_GRID" << endl;

    int ne = ma->GetNE();

    IntRange range = only_element >= 0 ? IntRange(only_element,only_element+1) : Range(ne);

    Array<Vec<D>> points;
    Array<Array<int>> cells;
    Array<int> cell_types;
    int pointcnt = 0;
    int total_cell_nums = 0;
    LocalHeap lh(1000000);
    for (int elnr : range)
    {
      HeapReset hr(lh);

      auto el = ma->GetElement(elnr);
      auto vertices = el.Vertices();

      Array<int> element_point_ids;
      for (int p : vertices)
      {
        points.Append(ma->GetPoint<D>(p));
        element_point_ids.Append(pointcnt++);
      }

      cells.Append(element_point_ids);
      cell_types.Append(ElementTypeToVTKType(el.GetType()));
    }

    ss << "POINTS " << points.Size() << " float" << endl;
    for (const Vec<D> & p : points)
    {
      for (int i = 0; i < D; ++i)
        ss << p[i] << " ";
      if (D==2)
        ss << "\t 0.0";
      ss << endl;
    }

    ss << "CELLS " << cells.Size() << " " << pointcnt + cells.Size() << endl;
    for (const Array<int> & c : cells)
    {
      ss << c.Size() << " ";
      for (int point : c)
        ss << point << " ";
      ss << endl;
    }

    ss << "CELL_TYPES " << cell_types.Size() << endl;
    for (int ct : cell_types)
      ss << ct << endl;

    ss << "CELL_DATA " << cells.Size() << endl;
    ss << "POINT_DATA " << points.Size() << endl;

    grid_str = ss.str();
  }
  /// Fill principil lattices (points and connections on subdivided reference simplex) in 2D
  template<>
  void VTKOutputQuad<2>::FillReferenceData(Array<IntegrationPoint> & ref_coords, Array<INT<2+1>> & ref_trigs)
  {
    enum { D = 2 };
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(1.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,1.0,0.0));
      ref_trigs.Append(INT<D+1>(0,1,2));
    }
    else
    {
      const int r = 1<<subdivision;
      const int s = r + 1;

      const double h = 1.0/r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i+j <= r; ++j)
          {
            ref_coords.Append(IntegrationPoint(j*h,i*h));
          }

      pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i+j <= r; ++j, pidx++)
          {
            // int pidx_curr = pidx;
            if (i+j == r) continue;
            int pidx_incr_i = pidx+1;
            int pidx_incr_j = pidx+s-i;

            ref_trigs.Append(INT<3>(pidx,pidx_incr_i,pidx_incr_j));

            int pidx_incr_ij = pidx_incr_j + 1;

            if(i+j+1<r)
              ref_trigs.Append(INT<3>(pidx_incr_i,pidx_incr_ij,pidx_incr_j));
          }
    }
  }

  /// Fill principil lattices (points and connections on subdivided reference simplex) in 3D
  template <int D>
  void VTKOutputQuad<D>::FillReferenceData(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_tets)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(1.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,1.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,0.0,1.0));
      ref_tets.Append(INT<D+1>(0,1,2,3));
    }
    else
    {
      const int r = 1<<subdivision;
      const int s = r + 1;

      const double h = 1.0/r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i+j <= r; ++j)
          for (int k = 0; i+j+k <= r; ++k)
          {
            ref_coords.Append(IntegrationPoint(i*h,j*h,k*h));
          }

      pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i+j <= r; ++j)
          for (int k = 0; i+j+k <= r; ++k, pidx++)
          {
            if (i+j+k == r) continue;
            // int pidx_curr = pidx;
            int pidx_incr_k = pidx+1;
            int pidx_incr_j = pidx+s-i-j;
            int pidx_incr_i = pidx+(s-i)*(s+1-i)/2-j;

            int pidx_incr_kj = pidx_incr_j + 1;

            int pidx_incr_ij = pidx+(s-i)*(s+1-i)/2-j + s-(i+1)-j;
            int pidx_incr_ki = pidx+(s-i)*(s+1-i)/2-j + 1;
            int pidx_incr_kij = pidx+(s-i)*(s+1-i)/2-j + s-(i+1)-j + 1;

            ref_tets.Append(INT<4>(pidx,pidx_incr_k,pidx_incr_j,pidx_incr_i));
            if (i+j+k+1 == r)
              continue;

            ref_tets.Append(INT<4>(pidx_incr_k,pidx_incr_kj,pidx_incr_j,pidx_incr_i));
            ref_tets.Append(INT<4>(pidx_incr_k,pidx_incr_kj,pidx_incr_ki,pidx_incr_i));

            ref_tets.Append(INT<4>(pidx_incr_j,pidx_incr_i,pidx_incr_kj,pidx_incr_ij));
            ref_tets.Append(INT<4>(pidx_incr_i,pidx_incr_kj,pidx_incr_ij,pidx_incr_ki));
           
            if (i+j+k+2 != r)
              ref_tets.Append(INT<4>(pidx_incr_kj,pidx_incr_ij,pidx_incr_ki,pidx_incr_kij));
          }              
    }
  }

  /// output of field data (coefficient values)
  template <int D> 
  void VTKOutputQuad<D>::PrintFieldData()
  {
    *fileout << "FIELD FieldData " << value_field.Size() << endl;

    for (auto field : value_field)
    {
      *fileout << field->Name() << " "
               << field->Dimension() << " "
               << field->Size()/field->Dimension() << " float" << endl;
      for (auto v : *field)
        *fileout << v << " ";
      *fileout << endl;
    }
  }

  template <int D>
  int VTKOutputQuad<D>::ElementTypeToVTKType(int et)
  {
    switch (et)
    {
    case ET_POINT: return VTK_VERTEX;
    case ET_SEGM: return VTK_LINE;
    case ET_TRIG: return VTK_TRIANGLE;
    case ET_QUAD: return VTK_QUAD;
    case ET_TET: return VTK_TETRA;
      // guessing here
    case ET_PRISM: return VTK_WEDGE;
    case ET_PYRAMID: return VTK_PYRAMID;
    case ET_HEX: return VTK_HEXAHEDRON;
    default: return -1;
    }
  }

  template <int D>
  void VTKOutputQuad<D>::Do (LocalHeap & lh, const BitArray * drawelems)
  {
    ostringstream filenamefinal;
    filenamefinal << filename;
    if (output_cnt > 0)
      filenamefinal << "_" << output_cnt;
    filenamefinal << ".vtk";
    fileout = make_shared<ofstream>(filenamefinal.str());
    cout << " Writing VTK-Output";
    if (output_cnt > 0)
      cout << " ( " << output_cnt << " )";
    cout << ":" << flush;

    output_cnt++;
    *fileout << grid_str;

    int ne = ma->GetNE();

    IntRange range = only_element >= 0 ? IntRange(only_element,only_element+1) : Range(ne);

    for ( int elnr : range)
    {
      if (drawelems && !(drawelems->Test(elnr)))
          continue;

      HeapReset hr(lh);

      ElementTransformation & eltrans = ma->GetTrafo (elnr, 0, lh);
      auto el = ma->GetElement(elnr);
      ELEMENT_TYPE et = el.GetType();
      const POINT3D *vertices = ElementTopology::GetVertices(et);
      int nv = ElementTopology::GetNVertices(et);

      for (int j = 0; j < nv; ++j)
      {
        MappedIntegrationPoint<D,D> mip(IntegrationPoint(vertices[j][0], vertices[j][1], vertices[j][2]), eltrans);
        for (int i = 0; i < coefs.Size(); ++i)
        {
          const int dim = coefs[i]->Dimension();
          FlatVector<> tmp(dim,lh);
          coefs[i]->Evaluate(mip,tmp);
          for (int d = 0; d < dim; ++d)
            value_field[i]->Append(tmp(d));
        }
      }

    }

    PrintFieldData();

    for (auto field : value_field)
      field->SetSize(0);

    cout << " Done." << endl;
  }

  // NumProcVTKOutput::NumProcVTKOutput (shared_ptr<PDE> apde, const Flags & flags)
  //   : NumProc (apde)
  // {
  //   const Array<string> & coefs_strings = flags.GetStringListFlag ("coefficients");

  //   Array<shared_ptr<CoefficientFunction>> coefs;
  //   for (int i = 0; i < coefs_strings.Size(); ++i)
  //     coefs.Append(apde->GetCoefficientFunction (coefs_strings[i]));

  //   if (apde->GetMeshAccess()->GetDimension() == 2)
  //     vtkout = make_shared<VTKOutput<2>>(coefs, flags, apde->GetMeshAccess());
  //   else
  //     vtkout = make_shared<VTKOutput<3>>(coefs, flags, apde->GetMeshAccess());
  // }


  // void NumProcVTKOutput::Do (LocalHeap & lh)
  // {
  //   vtkout->Do(lh);
  // }

  
  // static RegisterNumProc<NumProcVTKOutput> npvtkout("vtkoutput");

  template class VTKOutputQuad<2>;
  template class VTKOutputQuad<3>;
}

