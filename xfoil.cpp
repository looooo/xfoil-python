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

std::map<std::string, double> get_aero_components(const XFoil & foil)
{
    std::map<std::string, double> aero_coeff;
    aero_coeff["alfa"] = foil.alfa;
    aero_coeff["cl"] = foil.cl;
    aero_coeff["cm"] = foil.cm;
    aero_coeff["cd"] = foil.cd;
    aero_coeff["cdp"] = foil.cdp;
    aero_coeff["cdf"] = foil.cdf;
    aero_coeff["xcp"] = foil.xcp;
    return aero_coeff;
}

bool run(XFoil & foil, int max_iterations, bool initBL)
{
    foil.lwake = false;
    foil.lvconv = false;
    if(!foil.viscal())
    {
        foil.lvconv = false;
//      QString str =QObject::tr("CpCalc: local speed too large\n Compressibility corrections invalid");
        return false;
    }
    int i=0;
    while(i<max_iterations && !foil.lvconv)
    {
        if(foil.ViscousIter())
        {
            // if(m_x0 && m_y0)
            // {
            //  m_x0->append((double)i);
            //  m_y0->append(foil.rmsbl);
            // }
            // if(m_x1 && m_y1)
            // {
            //  m_x1->append((double)i);
            //  m_y1->append(foil.rmxbl);
            // }
            i++;
        }
        else i = max_iterations;
    }

    if(!foil.ViscalEnd())
    {
        foil.lvconv = false; //point is unconverged

        foil.setBLInitialized(false);
        foil.lipan = false;
        return true;// to exit loop
    }

    if(i>=max_iterations && !foil.lvconv)
    {
        if(initBL)
        {
            foil.setBLInitialized(false);
            foil.lipan = false;
        }
        foil.fcpmin();// Is it of any use ?
        return true;
    }
    if(!foil.lvconv)
    {
        foil.fcpmin();// Is it of any use ?
        return false;
    }
    else
    {
        //converged at last
        foil.fcpmin();// Is it of any use ?
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

Eigen::Matrix<double, Eigen::Dynamic, 2> get_normals(const XFoil& foil)
{
    Eigen::Matrix<double, Eigen::Dynamic, 2> normals;
    normals.resize(foil.nb, 2);
    for (int i = 0; i < foil.nb; i++)
    {
        normals(i, 0) = foil.nx[i + 1];
        normals(i, 1) = foil.ny[i + 1];
    }
    return normals;
}


Eigen::MatrixXd get_panel_data(const XFoil& foil)
{
    Eigen::MatrixXd panel_data;
    panel_data.resize(foil.nb, 4);
    for (int i = 0; i < foil.nb; i++)
    {
        panel_data(i, 0) = foil.cpi[i + 1];
        panel_data(i, 1) = foil.cpv[i + 1];
        panel_data(i, 2) = foil.qinv[i + 1];
        panel_data(i, 3) = foil.qvis[i + 1];
    }
    return panel_data;
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
            py::arg("bViscous")=true)
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
        .def("set_alfa", &XFoil::setAlpha)
        .def("set_cl", &XFoil::setClSpec)
        .def("specal", &XFoil::specal)
        .def("speccl", &XFoil::speccl)
        .def_readonly("numpoints", &XFoil::n)
        .def("get_coordiantes", &get_coordinates)
        .def("get_normals", &get_normals)
        .def("get_panel_data", &get_panel_data)
        .def_readwrite("lalfa", &XFoil::lalfa);

};
