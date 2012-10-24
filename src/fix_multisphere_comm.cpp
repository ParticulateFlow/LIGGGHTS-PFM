/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "fix_multisphere.h"
#include "comm.h"

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixMultisphere::pack_exchange(int i, double *buf)
{
    buf[0] = static_cast<double>(body_[i]);
    buf[1] = displace_[i][0];
    buf[2] = displace_[i][1];
    buf[3] = displace_[i][2];
    return 4;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixMultisphere::unpack_exchange(int nlocal, double *buf)
{
    body_[nlocal] = static_cast<int> (buf[0]);
    displace_[nlocal][0] = buf[1];
    displace_[nlocal][1] = buf[2];
    displace_[nlocal][2] = buf[3];
    return 4;
}

/* ----------------------------------------------------------------------
   forward comm
------------------------------------------------------------------------- */

void FixMultisphere::forward_comm()
{
    comm->forward_comm_fix(this);

    // force set of flag for next call
    fw_comm_flag_ = MS_COMM_UNDEFINED;
}

/* ----------------------------------------------------------------------
   pack comm
------------------------------------------------------------------------- */

int FixMultisphere::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
    /*NL*/ //fprintf(screen,"fw_comm_flag_ %d\n",fw_comm_flag_);

    if     (fw_comm_flag_ == MS_COMM_FW_BODY)
        return pack_comm_body(n,list,buf,pbc_flag,pbc);
    else if(fw_comm_flag_ == MS_COMM_FW_IMAGE_DISPLACE)
        return pack_comm_image_displace(n,list,buf,pbc_flag,pbc);
    else if(fw_comm_flag_ == MS_COMM_FW_V_OMEGA)
        return pack_comm_v_omega(n,list,buf,pbc_flag,pbc);
    else if(fw_comm_flag_ == MS_COMM_FW_F_TORQUE)
        return pack_comm_f_torque(n,list,buf,pbc_flag,pbc);
    else error->fix_error(FLERR,this,"FixMultisphere::pack_comm internal error");
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::pack_comm_body(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
    //we dont need to account for pbc here
    int i,j, m = 0;
    for (i = 0; i < n; i++)
    {
        j = list[i];

        buf[m++] = static_cast<double>(body_[j]);
    }
    return 1;
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::pack_comm_image_displace(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
    int *aimage = atom->image;

    //we dont need to account for pbc here
    int i,j, m = 0;
    for (i = 0; i < n; i++)
    {
        j = list[i];

        buf[m++] = static_cast<double>(aimage[j]);
        vectorToBuf3D(displace_[j],buf,m);
    }
    return 4;
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::pack_comm_v_omega(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
    double **v = atom->v;
    double **omega = atom->omega;

    //we dont need to account for pbc here
    int i,j, m = 0;
    for (i = 0; i < n; i++)
    {
        j = list[i];

        vectorToBuf3D(v[j],buf,m);
        vectorToBuf3D(omega[j],buf,m);
    }
    return 6;
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::pack_comm_f_torque(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
    double **f = atom->f;
    double **torque = atom->torque;

    //we dont need to account for pbc here
    int i,j, m = 0;
    int tag,flag;
    for (i = 0; i < n; i++)
    {
        j = list[i];

        //NP flag = 0 if single particle or body owned
        //NP flag = 1 if body NOT owned
        tag = body_[j];
        if(tag < 0) flag = 0;
        else flag = multisphere_.map(tag) < 0;
        /*NL*/// if(flag) fprintf(screen,"pack_comm atom %d with flag %d\n",atom->tag[j],flag);
        buf[m++] = static_cast<double>(flag);
        vectorToBuf3D(f[j],buf,m);
        vectorToBuf3D(torque[j],buf,m);
    }
    return 7;
}

/* ----------------------------------------------------------------------
   unpack comm
------------------------------------------------------------------------- */

void FixMultisphere::unpack_comm(int n, int first, double *buf)
{
    if     (fw_comm_flag_ == MS_COMM_FW_BODY)
        unpack_comm_body(n,first,buf);
    else if(fw_comm_flag_ == MS_COMM_FW_IMAGE_DISPLACE)
        unpack_comm_image_displace(n,first,buf);
    else if(fw_comm_flag_ == MS_COMM_FW_V_OMEGA)
        unpack_comm_v_omega(n,first,buf);
    else if(fw_comm_flag_ == MS_COMM_FW_F_TORQUE)
        unpack_comm_f_torque(n,first,buf);
    else error->fix_error(FLERR,this,"FixMultisphere::unpack_comm internal error");
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::unpack_comm_body(int n, int first, double *buf)
{
    int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++)
    {
        body_[i] = static_cast<int>(buf[m++]);
        /*NL*/ //fprintf(screen,"step %d: atom tag %d has body %d\n",update->ntimestep,atom->tag[i],body_[i]);
    }
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::unpack_comm_image_displace(int n, int first, double *buf)
{
    int i,m,last;
    int *aimage = atom->image;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++)
    {
        aimage[i] = static_cast<int>(buf[m++]);
        bufToVector3D(displace_[i],buf,m);
    }
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::unpack_comm_v_omega(int n, int first, double *buf)
{
    double **v = atom->v;
    double **omega = atom->omega;

    int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++)
    {
        bufToVector3D(v[i],buf,m);
        bufToVector3D(omega[i],buf,m);
    }
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::unpack_comm_f_torque(int n, int first, double *buf)
{
    int i,m,last,flag;
    double **f = atom->f;
    double **torque = atom->torque;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++)
    {
        //NP only apply force and torque if body NOT owned
        flag = static_cast<int>(buf[m++]);
        if(flag)
        {
            bufToVector3D(f[i],buf,m);
            bufToVector3D(torque[i],buf,m);
        }
        else m += 6;
    }
}

/* ----------------------------------------------------------------------
   forward comm
------------------------------------------------------------------------- */

void FixMultisphere::reverse_comm()
{
    comm->reverse_comm_fix(this);

    // force set of flag for next call
    rev_comm_flag_ = MS_COMM_UNDEFINED;
}

/* ----------------------------------------------------------------------
   pack reverse comm
------------------------------------------------------------------------- */

int FixMultisphere::pack_reverse_comm(int n, int first, double *buf)
{
    /*NL*/ //fprintf(screen,"rev_comm_flag_ %d\n",rev_comm_flag_);
    if     (rev_comm_flag_ == MS_COMM_REV_X_V_OMEGA)
        return pack_reverse_comm_x_v_omega(n,first,buf);
    else if(rev_comm_flag_ == MS_COMM_REV_V_OMEGA)
        return pack_reverse_comm_v_omega(n,first,buf);
    else error->fix_error(FLERR,this,"FixMultisphere::pack_reverse_comm internal error");
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::pack_reverse_comm_x_v_omega(int n, int first, double *buf)
{
    int i,m,last,tag,flag;

    double **x = atom->x;
    double **v = atom->v;
    double **omega = atom->omega;
    double *corner_ghost = fix_corner_ghost_->vector_atom;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {

        //NP flag = 0 if single particle or body NOT owned
        //NP flag = 1 if body owned
        tag = body_[i];

        if(tag < 0) flag = 0;
        else if(multisphere_.map(tag) >= 0) flag = 1;
        else if(corner_ghost[i] == 1.) flag = 1;
        else flag = 0;

        buf[m++] = static_cast<double>(flag);
        vectorToBuf3D(x[i],buf,m);
        vectorToBuf3D(v[i],buf,m);
        vectorToBuf3D(omega[i],buf,m);
    }
    return 10;
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::pack_reverse_comm_v_omega(int n, int first, double *buf)
{
    int i,m,last,tag,flag;

    double **v = atom->v;
    double **omega = atom->omega;
    double *corner_ghost = fix_corner_ghost_->vector_atom;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {

        //NP flag = 0 if single particle or body NOT owned
        //NP flag = 1 if body owned
        tag = body_[i];

        if(tag < 0) flag = 0;
        else if(multisphere_.map(tag) >= 0) flag = 1;
        else if(corner_ghost[i] == 1.) flag = 1;
        else flag = 0;

        buf[m++] = static_cast<double>(flag);
        vectorToBuf3D(v[i],buf,m);
        vectorToBuf3D(omega[i],buf,m);
    }
    return 7;
}

/* ----------------------------------------------------------------------
   unpack reverse comm
------------------------------------------------------------------------- */

void FixMultisphere::unpack_reverse_comm(int n, int *list, double *buf)
{
    if     (rev_comm_flag_ == MS_COMM_REV_X_V_OMEGA)
        unpack_reverse_comm_x_v_omega(n,list,buf);
    else if(rev_comm_flag_ == MS_COMM_REV_V_OMEGA)
        unpack_reverse_comm_v_omega(n,list,buf);
    else error->fix_error(FLERR,this,"FixMultisphere::unpack_reverse_comm internal error");
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::unpack_reverse_comm_x_v_omega(int n, int *list, double *buf)
{
    int i,j,flag,m = 0;

    int nlocal = atom->nlocal;
    double **x = atom->x;
    double **v = atom->v;
    double **omega = atom->omega;
    double *corner_ghost = fix_corner_ghost_->vector_atom;

    for (i = 0; i < n; i++) {
        j = list[i];

        flag = static_cast<int>(buf[m++]);
        if(flag)
        {
            bufToVector3D(x[j],buf,m);
            bufToVector3D(v[j],buf,m);
            bufToVector3D(omega[j],buf,m);
            //NP if unpacking at ghost, its a 'corner ghost', need to ensure that
            //NP information is transported across corners
            if(j >= nlocal)
                corner_ghost[j] = 1.;
        }
        else m += 9;
    }
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::unpack_reverse_comm_v_omega(int n, int *list, double *buf)
{
    int i,j,flag,m = 0;

    int nlocal = atom->nlocal;
    double **v = atom->v;
    double **omega = atom->omega;
    double *corner_ghost = fix_corner_ghost_->vector_atom;

    for (i = 0; i < n; i++) {
        j = list[i];

        flag = static_cast<int>(buf[m++]);
        if(flag)
        {
            bufToVector3D(v[j],buf,m);
            bufToVector3D(omega[j],buf,m);
            //NP if unpacking at ghost, its a 'corner ghost', need to ensure that
            //NP information is transported across corners
            if(j >= nlocal)
                corner_ghost[j] = 1.;
        }
        else m += 6;
    }
}

/* ----------------------------------------------------------------------
   restart functions
------------------------------------------------------------------------- */
//NP modified C.K.

int FixMultisphere::pack_restart(int, double *) { error->all(FLERR,"restart not implemented yet for multi-sphere feature"); return 0;}
void FixMultisphere::unpack_restart(int, int) {error->all(FLERR,"restart not implemented yet for multi-sphere feature");}
int FixMultisphere::size_restart(int) {error->all(FLERR,"restart not implemented yet for multi-sphere feature");return 0;}
int FixMultisphere::maxsize_restart() {error->all(FLERR,"restart not implemented yet for multi-sphere feature");return 0;}
