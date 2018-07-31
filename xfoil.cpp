#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include "XFoil.h"
#include <map>
#include <string>


namespace py = pybind11;

std::map<std::string, double> get_aero_components(const XFoil & xfoil)
{
    std::map<std::string, double> aero_coeff;
    aero_coeff["cl"] = xfoil.cl;
    aero_coeff["cm"] = xfoil.cm;
    aero_coeff["cd"] = xfoil.cd;
    aero_coeff["xcp"] = xfoil.xcp;
    return aero_coeff;
}

bool run(XFoil & xfoil, int max_iterations, bool initBL)
{
    if(!xfoil.viscal())
    {
        xfoil.lvconv = false;
//      QString str =QObject::tr("CpCalc: local speed too large\n Compressibility corrections invalid");
        return false;
    }
    int i=0;
    while(i<max_iterations && !xfoil.lvconv)
    {
        if(xfoil.ViscousIter())
        {
            // if(m_x0 && m_y0)
            // {
            //  m_x0->append((double)i);
            //  m_y0->append(xfoil.rmsbl);
            // }
            // if(m_x1 && m_y1)
            // {
            //  m_x1->append((double)i);
            //  m_y1->append(xfoil.rmxbl);
            // }
            i++;
        }
        else i = max_iterations;
    }

    if(!xfoil.ViscalEnd())
    {
        xfoil.lvconv = false; //point is unconverged

        xfoil.setBLInitialized(false);
        xfoil.lipan = false;
        return true;// to exit loop
    }

    if(i>=max_iterations && !xfoil.lvconv)
    {
        if(initBL)
        {
            xfoil.setBLInitialized(false);
            xfoil.lipan = false;
        }
        xfoil.fcpmin();// Is it of any use ?
        return true;
    }
    if(!xfoil.lvconv)
    {
        xfoil.fcpmin();// Is it of any use ?
        return false;
    }
    else
    {
        //converged at last
        xfoil.fcpmin();// Is it of any use ?
        return true;
    }
    return false;
}

Eigen::Matrix<double, Eigen::Dynamic, 2> get_coordinates(const XFoil& foil)
{
    Eigen::Matrix<double, Eigen::Dynamic, 2> coordinates;
    coordinates.resize(foil.nb, 2);
    for (int i = 0; i < foil.nb; i++)
    {
        coordinates(i, 0) = foil.xb[i + 1];
        coordinates(i, 1) = foil.yb[i + 1];
    }
    return coordinates;
}


PYBIND11_MODULE(xfoil, m) {
    m.doc() = "python bindings to xfoil";
    py::class_<XFoil>(m, "XFoil")
        .def(py::init<>())
        .def("initialize", &XFoil::initialize)
        .def("init_geometry", &XFoil::_initXFoilGeometry, py::arg("foil_coordinates"))
        .def("init_analysis", &XFoil::initXFoilAnalysis,
            py::arg("Re"), py::arg("alpha"), py::arg("Mach"), py::arg("NCrit"),
            py::arg("XtrTop")=1., py::arg("XtrBot")=1.,
            py::arg("reType")=1, py::arg("maType")=1, 
            py::arg("bViscous")=false)
        .def("get_aero_components", &get_aero_components)
        .def("run", &run,
            py::arg("max_iterations")=100,
            py::arg("initBL")=true)
        .def("viscal", &XFoil::viscal)
        .def("ViscousIter", &XFoil::ViscousIter)
        .def("ViscalEnd", &XFoil::ViscalEnd)
        .def_readonly("rmsbl", &XFoil::rmsbl)
        .def_readonly("lvconv", &XFoil::lvconv)
        .def_property("qinf", &XFoil::QInf, &XFoil::setQInf)
        .def("specal", &XFoil::specal)
        .def("speccl", &XFoil::speccl)
        .def_readonly("numpoints", &XFoil::nb)
        .def("get_coordiantes", &get_coordinates);
};
