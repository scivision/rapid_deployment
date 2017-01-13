/*
 * PyLpm.c -- Numerical python interface code for the Lag Profile Matrix
 *   (LPM) computation library of Open Radar
 *
 *   Copyright (c) 1999-2002, University of Tromsø and Massachusetts
 *   Institute of Technology, All Rights Reserved
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License as
 *    published by the Free Software Foundation; either version 2 of the
 *    License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *    USA
 */

/* $Id: PyLpm.c 5723 2011-03-09 14:05:00Z brideout $ */


#include <Python.h>
#include <numpy/oldnumeric.h>

#include "NCO.h"
#include "complex_type.h"
#include "mixer.h"
#include "filter.h"
#include "lpm.h"

/*****************************************************************************
 * Various signal-processing functions
 */

static char PyDSP_mixer_docstring[] =
    "interface to do_mixz_nco()\nmixer(z, samp_frq, osc_frq, phase)";
PyObject *
PyDSP_mixer(PyObject *self, PyObject *args)
{
    int ndims = 1, dims[1];

    double samp_frq, osc_frq, phase;
    NCO osc_sin, osc_cos;

    PyObject      *y     = NULL;
    PyArrayObject *z_in  = NULL;
    PyArrayObject *z_out = NULL;

    if (!PyArg_ParseTuple(args, "Oddd", &y, &samp_frq, &osc_frq, &phase))
	return NULL;

    /* Create oscillators */
    if (   nco_init(&osc_sin, osc_frq, 1.0/samp_frq) < 0
	|| nco_init(&osc_cos, osc_frq, 1.0/samp_frq) < 0) {
	PyErr_SetString(PyExc_ValueError, "frequency too high for oscillator");
	return NULL;
    }

    osc_sin.setphase(&osc_sin, phase);
    osc_cos.setphase(&osc_cos, phase + M_PI_2);

    if ((z_in = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
    {
	PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
	return NULL;
    }

    dims[0] = z_in->dimensions[0];

    if ((z_out = (PyArrayObject *)
	    PyArray_FromDims(ndims, dims, PyArray_CDOUBLE)) == NULL) {
	Py_DECREF(z_in);
	return NULL;
    }

    do_mixz_nco((complex *)z_out->data, (complex *)z_in->data, dims[0],
	&osc_sin, &osc_cos);

    Py_DECREF(z_in);

    return PyArray_Return(z_out);
}

static char PyDSP_mixvector_docstring[] =
    "create a fixed frq vector, suitable for mixing\n"
    "mixvector(len, samp_frq, osc_frq, phase)";
PyObject *
PyDSP_mixvector(PyObject *self, PyObject *args)
{
    int i, ndims = 1, dims[1];
    double samp_frq, osc_frq, phase;
    complex *vec;

    NCO osc_sin, osc_cos;
    PyArrayObject *out = NULL;

    if (!PyArg_ParseTuple(args, "iddd", &dims[0], &samp_frq, &osc_frq, &phase))
	return NULL;

    if (   nco_init(&osc_sin, osc_frq, 1.0/samp_frq) < 0
	|| nco_init(&osc_cos, osc_frq, 1.0/samp_frq) < 0) {
	PyErr_SetString(PyExc_ValueError, "frequency too high for oscillator");
	return NULL;
    }

    osc_sin.setphase(&osc_sin, phase);
    osc_cos.setphase(&osc_cos, phase + M_PI_2);

    if ((out = (PyArrayObject *)
		PyArray_FromDims(ndims, dims, PyArray_CDOUBLE)) == NULL) {
	return NULL;
    }

    vec = (complex *)out->data;

    for (i = 0; i < dims[0]; i++) {
	vec[i].re = osc_cos.next(&osc_cos);
	vec[i].im = osc_sin.next(&osc_sin);
    }

    return PyArray_Return(out);
}

static char PyDSP_filter_docstring[] =
    "interface to do_filtz()\nfilter(z, b, decimation)";
PyObject *
PyDSP_filter(PyObject *self, PyObject *args)
{
    int ndims = 1, dims[1];

    int decimation;
    complex *res;
    filt_t *f;

    PyObject      *x      = NULL;
    PyObject      *y      = NULL;
    PyArrayObject *z_in   = NULL;
    PyArrayObject *z_out  = NULL;
    PyArrayObject *filt_b = NULL;

    if (!PyArg_ParseTuple(args, "OOi", &x, &y, &decimation))
	return NULL;

    if ((z_in = (PyArrayObject *)
		PyArray_ContiguousFromObject(x, PyArray_CDOUBLE, 1, 1)) == NULL) {
	PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
	return NULL;
    }

    if ((filt_b = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_DOUBLE, 1, 1)) == NULL) {
	PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
	goto TYPE_ERROR_A;
    }

    if ((f = filter_from_array((double *)filt_b->data, filt_b->dimensions[0]))
	    == NULL) {
	PyErr_SetString(PyExc_RuntimeError, "Unable to create filter from input");
	goto FILTER_ERROR;
    }

    if ((dims[0] = do_filtz(&res, f, decimation, (complex *)z_in->data,
	    z_in->dimensions[0])) < 0) {
	PyErr_SetString(PyExc_RuntimeError, "Unable to filter input");
	goto FILTER_ERROR;
    }

    if ((z_out = (PyArrayObject *) PyArray_FromDimsAndData(ndims, dims,
		    PyArray_CDOUBLE, (char *)res)) == NULL)
	goto RESULT_ERROR;

    /* hand ownership of result data over to Python */
    z_out->flags |= OWN_DATA;

    filter_del(f);
    Py_DECREF(filt_b);
    Py_DECREF(z_in);

    return PyArray_Return(z_out);

RESULT_ERROR:
    free(res);
    filter_del(f);
FILTER_ERROR:
    Py_DECREF(filt_b);
TYPE_ERROR_A:
    Py_DECREF(z_in);
    return NULL;
}

/*****************************************************************************
 * Lag profile matrix functions
 */

static char PyLpm_startindex_docstring[] =
    "startindex(nzero, frac, ilag, flag) where flag = frational lag";
PyObject *
PyLpm_startindex(PyObject *self, PyObject *args)
{
    int nzero, frac, ilag, flag;

    if (!PyArg_ParseTuple(args, "iiii", &nzero, &frac, &ilag, &flag))
	return NULL;

    return Py_BuildValue("i", lp_startindex(nzero, frac, ilag, flag));
}

static char PyLpm_rstartindex_docstring[] =
    "rstartindex(nzero, frac, ilag, flag) where flag = frational lag";
PyObject *
PyLpm_rstartindex(PyObject *self, PyObject *args)
{
    int nzero, frac, ilag, flag;

    if (!PyArg_ParseTuple(args, "iiii", &nzero, &frac, &ilag, &flag))
	return NULL;

    return Py_BuildValue("i", rlp_startindex(nzero, frac, ilag, flag));
}

static char PyLpm_lplength_docstring[] =
    "lplength(nzero, frac, ilag, flag) where flag = frational lag";
PyObject *
PyLpm_lplength(PyObject *self, PyObject *args)
{
    int nzero, frac, ilag, flag;

    if (!PyArg_ParseTuple(args, "iiii", &nzero, &frac, &ilag, &flag))
	return NULL;

    return Py_BuildValue("i", lp_length(nzero, frac, ilag, flag));
}

static char PyLpm_length_docstring[] =
    "length(nzero, frac, ilag)";
PyObject *
PyLpm_length(PyObject *self, PyObject *args)
{
    int nzero, frac, ilag;

    if (!PyArg_ParseTuple(args, "iii", &nzero, &frac, &ilag))
	return NULL;

    return Py_BuildValue("i", lpm_length(nzero, frac, ilag));
}

static char PyLpm_rlength_docstring[] =
    "rlength(nzero, frac, ilag)";
PyObject *
PyLpm_rlength(PyObject *self, PyObject *args)
{
    int nzero, frac, ilag;

    if (!PyArg_ParseTuple(args, "iii", &nzero, &frac, &ilag))
	return NULL;

    return Py_BuildValue("i", rlpm_length(nzero, frac, ilag));
}

static char PyLpm_acfmac_docstring[] =
    "acfmac(lpm, nzero, nlags, data)";
PyObject *
PyLpm_acfmac(PyObject *self, PyObject *args)
{
    int nzero, nlags;
    PyObject *y;
    PyArrayObject *lpm = NULL, *z = NULL;

    if (!PyArg_ParseTuple(args, "O!iiO", &PyArray_Type, &lpm,
		&nzero, &nlags, &y))
	return NULL;

    if (!PyArray_Check(lpm) || lpm->descr->type_num != PyArray_CDOUBLE) {
	PyErr_SetString(PyExc_TypeError,
		"First argument to acfmac must be Numeric array "
		"of type complex double");
	return NULL;
    }

    if ((z = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    // printf("nzero:       %d\n", nzero);
    // printf("nlags:       %d\n", nlags);
    // printf("data length: %d\n", PyArray_Size((PyObject *)z));
    // printf("LPM length:  %d\n", PyArray_Size((PyObject *)lpm));
    lpm_acfmac((complex *)lpm->data, nzero, nlags, (complex *)z->data);

    Py_DECREF(z);
    return Py_BuildValue("");

TYPE_ERROR:
    Py_XDECREF(lpm);
    Py_XDECREF(z);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;
}

static char PyLpm_decode_docstring[] =
    "decode(rlpm, code, weight, nzero, frac, nilags, lpm)";
PyObject *
PyLpm_decode(PyObject *self, PyObject *args)
{
    char *code;
    int nzero, nilags, frac;
    double weight;
    PyObject *y;
    PyArrayObject *rlpm = NULL, *lpm = NULL;

    if (!PyArg_ParseTuple(args, "O!sdiiiO", &PyArray_Type, &rlpm,
		&code, &weight, &nzero, &frac, &nilags, &y))
	return NULL;

    // printf("code: %s\n", code);

    if (!PyArray_Check(rlpm) || rlpm->descr->type_num != PyArray_CDOUBLE) {
	PyErr_SetString(PyExc_TypeError,
		"First argument to decode must be Numeric array of type float");
	return NULL;
    }

    if ((lpm = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    /*
    printf("PyLpm_decode: rlpm: %d (%p)  code: <%s>  nzero: %d  frac: %d"
	    "  nilags: %d\n\tlpm: %d (%p)\n",
	    PyArray_Size(rlpm), rlpm->data, code, nzero, frac,
	    nilags, PyArray_Size(lpm), lpm->data);
     */

    rlpm_decode((complex *)rlpm->data, code, weight, nzero, frac, nilags,
	    (complex *)lpm->data);

    Py_DECREF(lpm);
    return Py_BuildValue("");

TYPE_ERROR:
    Py_XDECREF(rlpm);
    Py_XDECREF(lpm);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;

}



static char PyLpm_decode_walsh_docstring[] =
    "decode_walsh(rlpm, code_value, code_id, weight, nzero, nilags, lpm)";
PyObject *
PyLpm_decode_walsh(PyObject *self, PyObject *args)
{
    int code_id, code_size;
    int nzero, nilags;
    double weight;
    PyObject *y;
    PyArrayObject *code_val = NULL;
    PyArrayObject *rlpm = NULL, *lpm = NULL;


    if (!PyArg_ParseTuple(args, "O!O!idiiO", &PyArray_Type, &rlpm,
		&PyArray_Type, &code_val, &code_id, &weight, &nzero, &nilags, &y))
	return NULL;


    if (!PyArray_Check(rlpm) || rlpm->descr->type_num != PyArray_CDOUBLE) {
	PyErr_SetString(PyExc_TypeError,
		"First argument to decode_walsh must be Numeric array of type float");
	return NULL;
    }

    if (!PyArray_Check(code_val) || code_val->descr->type_num != PyArray_INT) {
	PyErr_SetString(PyExc_TypeError,
		"Second argument to decode_walsh a must be Numeric array of type int");
	return NULL;
    }


    if ((lpm = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    /*
    printf("PyLpm_decode_walsh: rlpm: %d (%p)  code_val: %d (%p) code_id: <%s>  nzero: %d"
	    "  nilags: %d\n\tlpm: %d (%p)\n",
	    PyArray_Size(rlpm), rlpm->data, PyArray_Size(code_val), code_val->data, code_id, nzero,
	    nilags, PyArray_Size(lpm), lpm->data);
     */

    code_size = PyArray_Size((PyObject *)code_val);

    rlpm_decode_walsh((complex *)rlpm->data, (int *)code_val->data, code_size,
                      code_id, weight, nzero, nilags,(complex *)lpm->data);


    Py_DECREF(lpm);
    return Py_BuildValue("");

TYPE_ERROR:
    Py_XDECREF(rlpm);
    Py_XDECREF(lpm);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;

}

static char PyLpm_acfmatrix_trapezoidal_docstring[] =
    "acfmatrix_trapezoidal(lpm, nzero, frac, nilags, gating, overlap)";
PyObject *
PyLpm_acfmatrix_trapezoidal(PyObject *self, PyObject *args)
{
    PyArrayObject *lpm = NULL;
    PyArrayObject *acf = NULL;
    PyObject      *y;

    int nzero, frac, nilags, gating, overlap;
    int _ngates, ngates, dims[2];

    if (!PyArg_ParseTuple(args, "O!iiiii", &PyArray_Type, &y,
		&nzero, &frac, &nilags, &gating, &overlap))
	return NULL;

    if ((lpm = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    ngates = n_igates_trapezoidal(nzero, nilags, gating, overlap);

    dims[0] = ngates;
    dims[1] = nilags*frac;
    /* must match nlags in function lpm_acfmatrix_trapezoidal_na in file lpm.c */

    if ((acf = (PyArrayObject *)
		PyArray_FromDims(2, dims, PyArray_CDOUBLE)) == NULL)
	return NULL;

    /*
    printf("acf: %p lpm: %p nzero: %d frac: %d nilags: %d\n"
	    "ngates: %d gating: %d overlap: %d\n",
	    acf->data, lpm->data, nzero, frac, nilags, ngates, gating, overlap);
     */

    /* Populate the new 2D matrix with ACFs */
    _ngates = lpm_acfmatrix_trapezoidal_na((complex *)acf->data, (complex *)lpm->data,
					   nzero, frac, nilags, ngates, gating, overlap);

    Py_DECREF(lpm);
    return PyArray_Return(acf);

TYPE_ERROR:
    Py_XDECREF(lpm);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;
}


static char PyLpm_acfmatrix_inversetrapezoidal_docstring[] =
    "acfmatrix_inversetrapezoidal(lpm, nzero, frac, nilags, gating, overlap)";
PyObject *
PyLpm_acfmatrix_inversetrapezoidal(PyObject *self, PyObject *args)
{
    PyArrayObject *lpm = NULL;
    PyArrayObject *acf = NULL;
    PyObject      *y;

    int nzero, frac, nilags, gating, overlap;
    int _ngates, ngates, dims[2];

    if (!PyArg_ParseTuple(args, "O!iiiii", &PyArray_Type, &y,
		&nzero, &frac, &nilags, &gating, &overlap))
	return NULL;

    if ((lpm = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    ngates = n_igates_inversetrapezoidal(nzero, nilags, gating, overlap);

    dims[0] = ngates;
    dims[1] = nilags*frac;
    /* must match nlags in function lpm_acfmatrix_inversetrapezoidal_na in file lpm.c */

    if ((acf = (PyArrayObject *)
		PyArray_FromDims(2, dims, PyArray_CDOUBLE)) == NULL)
	return NULL;

    /*
    printf("acf: %p lpm: %p nzero: %d frac: %d nilags: %d\n"
	    "ngates: %d gating: %d overlap: %d\n",
	    acf->data, lpm->data, nzero, frac, nilags, ngates, gating, overlap);
     */

    /* Populate the new 2D matrix with ACFs */
    _ngates = lpm_acfmatrix_inversetrapezoidal_na((complex *)acf->data, (complex *)lpm->data,
						  nzero, frac, nilags, ngates, gating, overlap);

    Py_DECREF(lpm);
    return PyArray_Return(acf);

TYPE_ERROR:
    Py_XDECREF(lpm);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;
}

static char PyLpm_racfmatrix_docstring[] =
    "racfmatrix(rlpm, nzero, frac, nilags, gating, overlap)";
PyObject *
PyLpm_racfmatrix(PyObject *self, PyObject *args)
{
    PyArrayObject *rlpm = NULL;
    PyArrayObject *acf = NULL;
    PyObject      *y;

    int nzero, frac, nilags, gating, overlap;
    int _ngates, ngates, dims[2];

    if (!PyArg_ParseTuple(args, "O!iiiii", &PyArray_Type, &y,
		&nzero, &frac, &nilags, &gating, &overlap))
	return NULL;

    if ((rlpm = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    ngates = n_fgates(nzero, nilags, frac, gating, overlap);

    dims[0] = ngates;
    dims[1] = (nilags-1)*frac + 1;
    /* must match nlags in function rlpm_acfmatrix_na in file lpm.c */

    if ((acf = (PyArrayObject *)
		PyArray_FromDims(2, dims, PyArray_CDOUBLE)) == NULL)
	return NULL;

    /*
    printf("acf: %p rlpm: %p nzero: %d frac: %d nilags: %d\n"
	    "ngates: %d gating: %d overlap: %d\n",
	    acf->data, rlpm->data, nzero, frac, nilags, ngates, gating, overlap);
     */

    /* Populate the new 2D matrix with ACFs */
    _ngates = rlpm_acfmatrix_na((complex *)acf->data, (complex *)rlpm->data,
	nzero, frac, nilags, ngates, gating, overlap);

    Py_DECREF(rlpm);
    return PyArray_Return(acf);

TYPE_ERROR:
    Py_XDECREF(rlpm);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;
}


/*****************************************************************************
 * Cross-correlation Lag profile matrix functions
 */

static char PyLpm_xstartindex_docstring[] =
    "xstartindex(nzero, frac, nilags, ilag, flag) where flag = frational lag";
PyObject *
PyLpm_xstartindex(PyObject *self, PyObject *args)
{
    int nzero, frac, nilags, ilag, flag;

    if (!PyArg_ParseTuple(args, "iiiii", &nzero, &frac, &nilags, &ilag, &flag))
	return NULL;

    return Py_BuildValue("i", xlp_startindex(nzero, frac, nilags, ilag, flag));
}

static char PyLpm_rxstartindex_docstring[] =
    "rxstartindex(nzero, frac, nilags, ilag, flag) where flag = frational lag";
PyObject *
PyLpm_rxstartindex(PyObject *self, PyObject *args)
{
    int nzero, frac, nilags, ilag, flag;

    if (!PyArg_ParseTuple(args, "iiiii", &nzero, &frac, &nilags, &ilag, &flag))
	return NULL;

    return Py_BuildValue("i", rxlp_startindex(nzero, frac, nilags, ilag, flag));
}

static char PyLpm_xlength_docstring[] =
    "xlength(nzero, frac, nilags)";
PyObject *
PyLpm_xlength(PyObject *self, PyObject *args)
{
    int nzero, frac, nilags;

    if (!PyArg_ParseTuple(args, "iii", &nzero, &frac, &nilags))
	return NULL;

    return Py_BuildValue("i", xlpm_length(nzero, frac, nilags));
}

static char PyLpm_rxlength_docstring[] =
    "rxlength(nzero, frac, nilags)";
PyObject *
PyLpm_rxlength(PyObject *self, PyObject *args)
{
    int nzero, frac, nilags;

    if (!PyArg_ParseTuple(args, "iii", &nzero, &frac, &nilags))
	return NULL;

    return Py_BuildValue("i", rxlpm_length(nzero, frac, nilags));
}

static char PyLpm_xmac_docstring[] =
    "xmac(xlpm, nzero, nlags, data1, data2)";
PyObject *
PyLpm_xmac(PyObject *self, PyObject *args)
{
    int nzero, nlags;
    PyObject *y1, *y2;
    PyArrayObject *xlpm = NULL, *z1 = NULL, *z2 = NULL;

    if (!PyArg_ParseTuple(args, "O!iiOO", &PyArray_Type, &xlpm,
		&nzero, &nlags, &y1, &y2))
	return NULL;

    if (!PyArray_Check(xlpm) || xlpm->descr->type_num != PyArray_CDOUBLE) {
	PyErr_SetString(PyExc_TypeError,
		"First argument to xmac must be Numeric array of type "
		"complex double");
	return NULL;
    }

    if ((z1 = (PyArrayObject *)
	    PyArray_ContiguousFromObject(y1, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR1;

    if ((z2 = (PyArrayObject *)
	    PyArray_ContiguousFromObject(y2, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR2;

    // printf("nzero:       %d\n", nzero);
    // printf("nlags:       %d\n", nlags);
    // printf("data length: %d\n", PyArray_Size((PyObject *)z));
    // printf("LPM length:  %d\n", PyArray_Size((PyObject *)lpm));
    xlpm_mac((complex *)xlpm->data, nzero, nlags,
	    (complex *)z1->data, (complex *)z2->data);

    Py_DECREF(z1);
    Py_DECREF(z2);
    return Py_BuildValue("");

TYPE_ERROR2:
    Py_XDECREF(z2);
TYPE_ERROR1:
    Py_XDECREF(z1);
    Py_XDECREF(xlpm);

    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;
}

static char PyLpm_xdecode_docstring[] =
    "xdecode(rxlpm, code, weight, nzero, frac, nilags, xlpm)";
PyObject *
PyLpm_xdecode(PyObject *self, PyObject *args)
{
    char *code;
    int nzero, nilags, frac;
    double weight;
    PyObject *y;
    PyArrayObject *rxlpm = NULL, *xlpm = NULL;

    if (!PyArg_ParseTuple(args, "O!sdiiiO", &PyArray_Type, &rxlpm,
		&code, &weight, &nzero, &frac, &nilags, &y))
	return NULL;

    // printf("code: %s\n", code);

    if (!PyArray_Check(rxlpm) || rxlpm->descr->type_num != PyArray_CDOUBLE) {
	PyErr_SetString(PyExc_TypeError,
		"First argument to acfmac must be Numeric array of type float");
	return NULL;
    }

    if ((xlpm = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    /*
    printf("PyLpm_decode: rxlpm: %d (%p)  code: <%s>  nzero: %d  frac: %d"
	    "  nilags: %d\n\txlpm: %d (%p)\n",
	    PyArray_Size(rxlpm), rxlpm->data, code, nzero, frac, nilags,
	    PyArray_Size(xlpm),   xlpm->data);
     */

    rxlpm_decode((complex *)rxlpm->data, code, weight, nzero, frac, nilags,
	    (complex *)xlpm->data);

    Py_DECREF(xlpm);
    return Py_BuildValue("");

TYPE_ERROR:
    Py_XDECREF(rxlpm);
    Py_XDECREF(xlpm);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;

}

static char PyLpm_xcfmatrix_trapezoidal_docstring[] =
    "xcfmatrix_trapezoidal(xlpm, nzero, frac, nilags, gating, overlap)";
PyObject *
PyLpm_xcfmatrix_trapezoidal(PyObject *self, PyObject *args)
{
    PyArrayObject *xlpm = NULL;
    PyArrayObject *xcf = NULL;
    PyObject      *y;

    int nzero, frac, nilags, gating, overlap;
    int _ngates, ngates, dims[2];

    if (!PyArg_ParseTuple(args, "O!iiiii", &PyArray_Type, &y,
		&nzero, &frac, &nilags, &gating, &overlap))
	return NULL;

    if ((xlpm = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    ngates = n_igates_trapezoidal(nzero, nilags, gating, overlap);

    dims[0] = ngates;
    dims[1] = 2*nilags - 1;
    /* must match nlags in function xlpm_xcfmatrix_na in file lpm.c */

    if ((xcf = (PyArrayObject *)
		PyArray_FromDims(2, dims, PyArray_CDOUBLE)) == NULL)
	return NULL;

    /*
    printf("xcf: %p xlpm: %p nzero: %d frac: %d nilags: %d\n"
	    "ngates: %d gating: %d overlap: %d\n", xcf->data, xlpm->data,
	    nzero, frac, nilags, ngates, gating, overlap);
     */

    /* Populate the new 2D matrix with ACFs */
    _ngates = xlpm_xcfmatrix_trapezoidal_na((complex *)xcf->data, (complex *)xlpm->data,
	nzero, frac, nilags, ngates, gating, overlap);

    Py_DECREF(xlpm);
    return PyArray_Return(xcf);

TYPE_ERROR:
    Py_XDECREF(xlpm);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;
}

static char PyLpm_rxcfmatrix_docstring[] =
    "rxcfmatrix(rxlpm, nzero, frac, nilags, gating, overlap)";
PyObject *
PyLpm_rxcfmatrix(PyObject *self, PyObject *args)
{
    PyArrayObject *rxlpm = NULL;
    PyArrayObject *xcf   = NULL;
    PyObject      *y;

    int nzero, frac, nilags, gating, overlap;
    int _ngates, ngates, dims[2];

    if (!PyArg_ParseTuple(args, "O!iiiii", &PyArray_Type, &y,
		&nzero, &frac, &nilags, &gating, &overlap))
	return NULL;

    if ((rxlpm = (PyArrayObject *)
		PyArray_ContiguousFromObject(y, PyArray_CDOUBLE, 1, 1)) == NULL)
	goto TYPE_ERROR;

    ngates = n_fgates(nzero, nilags, frac, gating, overlap);

    dims[0] = ngates;
    dims[1] = 2*(nilags-1)*frac + 1;
    /* must match nlags in function rlpm_acfmatrix_na in file lpm.c */

    if ((xcf = (PyArrayObject *)
		PyArray_FromDims(2, dims, PyArray_CDOUBLE)) == NULL)
	return NULL;

    /*
    printf("acf: %p rlpm: %p nzero: %d frac: %d nilags: %d\n"
	    "ngates: %d gating: %d overlap: %d\n",
	    acf->data, rlpm->data, nzero, frac, nilags, ngates, gating, overlap);
     */

    /* Populate the new 2D matrix with ACFs */
    _ngates = rxlpm_xcfmatrix_na((complex *)xcf->data, (complex *)rxlpm->data,
	nzero, frac, nilags, ngates, gating, overlap);

    Py_DECREF(rxlpm);
    return PyArray_Return(xcf);

TYPE_ERROR:
    Py_XDECREF(rxlpm);
    PyErr_SetString(PyExc_TypeError, "Inappropriate type(s) for input(s)");
    return NULL;
}


void
initPyLpm(void)
{
    static PyMethodDef PyLpm_methods[] = {
	{ "mixer",		PyDSP_mixer,		METH_VARARGS,
	    PyDSP_mixer_docstring },
	{ "mixvector",		PyDSP_mixvector,	METH_VARARGS,
	    PyDSP_mixvector_docstring },
	{ "filter",		PyDSP_filter,		METH_VARARGS,
	    PyDSP_filter_docstring },
	{ "startindex",		PyLpm_startindex,	METH_VARARGS,
	    PyLpm_startindex_docstring },
	{ "rstartindex",	PyLpm_rstartindex,	METH_VARARGS,
	    PyLpm_rstartindex_docstring },
	{ "lplength",           PyLpm_lplength,         METH_VARARGS,
	    PyLpm_lplength_docstring },
	{ "length",		PyLpm_length,		METH_VARARGS,
	    PyLpm_length_docstring },
	{ "rlength",		PyLpm_rlength,		METH_VARARGS,
	    PyLpm_rlength_docstring },
	{ "acfmac",		PyLpm_acfmac,		METH_VARARGS,
	    PyLpm_acfmac_docstring },
	{ "decode",		PyLpm_decode,		METH_VARARGS,
	    PyLpm_decode_docstring },
	{ "decode_walsh",       PyLpm_decode_walsh,     METH_VARARGS,
	    PyLpm_decode_walsh_docstring },
	{ "acfmatrix_trapezoidal", PyLpm_acfmatrix_trapezoidal,	METH_VARARGS,
	    PyLpm_acfmatrix_trapezoidal_docstring },
	{ "acfmatrix_inversetrapezoidal", PyLpm_acfmatrix_inversetrapezoidal,	METH_VARARGS,
	    PyLpm_acfmatrix_inversetrapezoidal_docstring },
	{ "racfmatrix",		PyLpm_racfmatrix,	METH_VARARGS,
	    PyLpm_racfmatrix_docstring },
	{ "xstartindex",	PyLpm_xstartindex,	METH_VARARGS,
	    PyLpm_xstartindex_docstring },
	{ "rxstartindex",	PyLpm_rxstartindex,	METH_VARARGS,
	    PyLpm_rxstartindex_docstring },
	{ "xlength",		PyLpm_xlength,		METH_VARARGS,
	    PyLpm_xlength_docstring },
	{ "rxlength",		PyLpm_rxlength,		METH_VARARGS,
	    PyLpm_rxlength_docstring },
	{ "xmac",		PyLpm_xmac,		METH_VARARGS,
	    PyLpm_xmac_docstring },
	{ "xdecode",		PyLpm_xdecode,		METH_VARARGS,
	    PyLpm_xdecode_docstring },
	{ "xcfmatrix_trapezoidal", PyLpm_xcfmatrix_trapezoidal,	METH_VARARGS,
	    PyLpm_xcfmatrix_trapezoidal_docstring },
	{ "rxcfmatrix",		PyLpm_rxcfmatrix,	METH_VARARGS,
	    PyLpm_rxcfmatrix_docstring },
	{ NULL, NULL },
    };

    Py_InitModule("PyLpm", PyLpm_methods);

    // Without the following call, any call to Numeric functions will
    // crash Python!
    import_array();
}

