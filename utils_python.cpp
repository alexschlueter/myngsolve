#ifdef NGS_PYTHON

#include <python_ngstd.hpp>
#include "vtkoutputquad.hpp"
#include "randomcf.hpp"

using namespace ngfem;

void ExportNgsAppsUtils() 
{
  cout << "exporting ngsapps.utils"  << endl;

  // Export RandomCoefficientFunction to python (name "RandomCF")

  bp::class_<RandomCoefficientFunction, shared_ptr<RandomCoefficientFunction>, bp::bases<CoefficientFunction>, boost::noncopyable>("RandomCF", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](double lower, double upper)
                           {
                             return make_shared<RandomCoefficientFunction> (lower, upper);
                           }),
          bp::default_call_policies(),        // need it to use named arguments
          (bp::arg("lower_bound")=0.0,bp::arg("upper_bound")=1.0)
           ));
  
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<RandomCoefficientFunction>);

  bp::implicitly_convertible 
    <shared_ptr<RandomCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 

  using namespace ngcomp;
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<BaseVTKOutput>);
  bp::class_<BaseVTKOutput, shared_ptr<BaseVTKOutput>,  boost::noncopyable>("VTKOutputQuad", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](shared_ptr<MeshAccess> ma, bp::list coefs_list,
                              bp::list names_list, string filename, int subdivision, int only_element)
                           {
                             Array<shared_ptr<CoefficientFunction> > coefs
                               = makeCArray<shared_ptr<CoefficientFunction>> (coefs_list);
                             Array<string > names
                               = makeCArray<string> (names_list);
                             shared_ptr<BaseVTKOutput> ret;
                             if (ma->GetDimension() == 2)
                               ret = make_shared<VTKOutputQuad<2>> (ma, coefs, names, filename, subdivision, only_element);
                             else
                               ret = make_shared<VTKOutputQuad<3>> (ma, coefs, names, filename, subdivision, only_element);
                             return ret;
                           }),

          bp::default_call_policies(),     // need it to use named arguments
          (
            bp::arg("ma"),
            bp::arg("coefs")= bp::list(),
            bp::arg("names") = bp::list(),
            bp::arg("filename") = "vtkout",
            bp::arg("subdivision") = 0,
            bp::arg("only_element") = -1
            )
           )
      )

    .def("Do", FunctionPointer([](BaseVTKOutput & self, int heapsize)
                               { 
                                 LocalHeap lh (heapsize, "VTKOutput-heap");
                                 self.Do(lh);
                               }),
         (bp::arg("self"),bp::arg("heapsize")=1000000))
    .def("Do", FunctionPointer([](BaseVTKOutput & self, const BitArray * drawelems, int heapsize)
                               { 
                                 LocalHeap lh (heapsize, "VTKOutput-heap");
                                 self.Do(lh, drawelems);
                               }),
         (bp::arg("self"),bp::arg("drawelems"),bp::arg("heapsize")=1000000))
    
    ;
  
}

BOOST_PYTHON_MODULE(libngsapps_utils) 
{
  ExportNgsAppsUtils();
}
#endif // NGSX_PYTHON
