//                                Finite Element Assignment 4
//
// Description:  You will implement the shape functions and the derivatives in this header file.
//               This header file will then be used in subsequent .cpp files as part of your assignment.
//               You will need to add your code in the sections marked TODO.
//               At the minimum you will need to implement
//               a) shape functions/their derivative for Segments in
//                  "LG_SegmentElement::CalcShape" and "LG_SegmentElement::CalcDShape"
//               b) shape functions/their derivative for Quads in
//                  "LG_QuadrilateralElement::CalcShape" and "LG_QuadrilateralElement::CalcDShape"
//               c) THIS IS NOT PART OF THE ASSIGNMENT but all the infrastructure is in place also for Triangles.
//                  So if you like, you can implement shape functions/their derivatives for Triangles in
//                  "LG_TriangleElement::CalcShape" and "LG_TriangleElement::CalcDShape"
//                  That will enable you to run the rest of this assignment on triangular meshes as well.
//

#ifndef LAGRANGE_ELEMENTS
#define LAGRANGE_ELEMENTS

#include <mfem.hpp>
#include <fstream>
#include <iostream>
#include <queue>

using namespace std;
using namespace mfem;

class LG_SegmentElement : public NodalTensorFiniteElement
{
public:
  LG_SegmentElement(
      const int p,
      const int btype = BasisType::GaussLobatto);
  virtual void CalcShape(
      const IntegrationPoint &ip,
      Vector &shape) const;
  virtual void CalcDShape(
      const IntegrationPoint &ip,
      DenseMatrix &dshape) const;
};

class LG_QuadrilateralElement : public NodalTensorFiniteElement
{
public:
  LG_QuadrilateralElement(
      const int p,
      const int btype = BasisType::GaussLobatto);
  virtual void CalcShape(
      const IntegrationPoint &ip,
      Vector &shape) const;
  virtual void CalcDShape(
      const IntegrationPoint &ip,
      DenseMatrix &dshape) const;
};


class LG_TriangleElement : public NodalFiniteElement
{
public:
  LG_TriangleElement(
      const int p,
      const int btype = BasisType::GaussLobatto);
  virtual void CalcShape(
      const IntegrationPoint &ip,
      Vector &shape) const;
  virtual void CalcDShape(
      const IntegrationPoint &ip,
      DenseMatrix &dshape) const;
};

class LG_FECollection : public FiniteElementCollection
{
protected:
  int b_type;
  char lg_name[32];
  FiniteElement *LG_Elements[Geometry::NumGeom];
  int LG_dof[Geometry::NumGeom];
  int *SegDofOrd[2], *TriDofOrd[6], *QuadDofOrd[8], *TetDofOrd[24];

public:
  explicit LG_FECollection(
      const int p,
      const int dim = 3,
      const int btype = BasisType::GaussLobatto);

  virtual const FiniteElement *FiniteElementForGeometry(
     Geometry::Type GeomType) const
  { return LG_Elements[GeomType]; }
  virtual int DofForGeometry(Geometry::Type GeomType) const
  { return LG_dof[GeomType]; }
  virtual const int *DofOrderForOrientation(
      Geometry::Type GeomType,
      int Or) const;
  virtual const char *Name() const
  { return lg_name; }
  virtual int GetContType() const
  { return CONTINUOUS; }
  FiniteElementCollection *GetTraceCollection() const;

  int GetBasisType() const
  { return b_type; }
  /// Get the Cartesian to local LG dof map
  const int *GetDofMap(Geometry::Type GeomType) const;

  virtual ~LG_FECollection();
};


// LG_SegmentElement implementation
LG_SegmentElement::LG_SegmentElement(
    const int p,
    const int btype)
   : NodalTensorFiniteElement(1, p, VerifyClosed(btype), H1_DOF_MAP)
{
  const double *cp = poly1d.ClosedPoints(p, b_type);

  Nodes.IntPoint(0).x = cp[0];
  Nodes.IntPoint(1).x = cp[p];
  for (int i = 1; i < p; i++)
  {
    Nodes.IntPoint(i+1).x = cp[i];
  }
}

void LG_SegmentElement::CalcShape(
    const IntegrationPoint &ip,
    Vector &shape) const
{
  // get the order from the class member "order"
  const int p = order;
  // Note: in all cases boundary nodes are the first two and then the mid nodes
  if (p == 1)
  {
    // TODO: implement 1st order shapes for segments
    // Note the parametric coordinate is ip.x
    shape(0) = 0.;
    shape(1) = 0.;
  }
  else if (p == 2)
  {
    // TODO: implement 2nd order shapes for segments
    shape(0) = 0.;
    shape(1) = 0.;
    shape(2) = 0.;
  }
  else
  {
    MFEM_VERIFY(0, "unimplemented order");
  }
}

void LG_SegmentElement::CalcDShape(
    const IntegrationPoint &ip,
    DenseMatrix &dshape) const
{
  // get the order from the class member "order"
  const int p = order;
  // Note: in all cases boundary nodes are the first two and then the mid nodes
  if (p == 1)
  {
    // TODO: implement first derivative of 1st order shapes for segments
    dshape(0,0) = 0.;
    dshape(1,0) = 0.;
  }
  else if (p == 2)
  {
    // TODO: implement first derivative of 2st order shapes for segments
    dshape(0,0) = 0.;
    dshape(1,0) = 0.;
    dshape(2,0) = 0.;
  }
  else
  {
    MFEM_VERIFY(0, "unimplemented order");
  }
}

// LG_QuadrilateralElement implementation
LG_QuadrilateralElement::LG_QuadrilateralElement(
    const int p,
    const int btype)
   : NodalTensorFiniteElement(2, p, VerifyClosed(btype), H1_DOF_MAP)
{
  const double *cp = poly1d.ClosedPoints(p, b_type);

  int o = 0;
  for (int j = 0; j <= p; j++)
    for (int i = 0; i <= p; i++)
      Nodes.IntPoint(dof_map[o++]).Set2(cp[i], cp[j]);
}

void LG_QuadrilateralElement::CalcShape(
    const IntegrationPoint &ip,
    Vector &shape) const
{
  const int p = order;

  // Note: In MFEM quads are defined on [0,1]x[0,1]
  //       Nodes for 1st order quads are as follows
  //       node 0: (0,0)
  //       node 1: (1,0)
  //       node 2: (1,1)
  //       node 3: (0,1)
  //
  //       Nodes for 2nd order quads are as follows
  //       node 0: (0,0)
  //       node 1: (1,0)
  //       node 2: (1,1)
  //       node 3: (0,1)
  //       node 4: (1/2,0)
  //       node 5: (1,1/2)
  //       node 6: (1/2,1)
  //       node 7: (0,1/2)
  //       node 8: (1/2,1/2)
  if (p == 1)
  {
    // TODO: implement 1st order shapes for quads
    shape(0) = 0.;
    shape(1) = 0.;
    shape(2) = 0.;
    shape(3) = 0.;
  }
  else if (p == 2)
  {
    // TODO: implement 2st order shapes for quads (9-nodes)
    shape(0) = 0.;
    shape(1) = 0.;
    shape(2) = 0.;
    shape(3) = 0.;
    shape(4) = 0.;
    shape(5) = 0.;
    shape(6) = 0.;
    shape(7) = 0.;
    shape(8) = 0.;
  }
  else
  {
    MFEM_VERIFY(0, "unimplemented order");
  }

}

void LG_QuadrilateralElement::CalcDShape(
    const IntegrationPoint &ip,
    DenseMatrix &dshape) const
{
  const int p = order;

  // Note: In MFEM quads are defined on [0,1]x[0,1]
  //       Nodes for 1st order quads are as follows
  //       node 0: (0,0)
  //       node 1: (1,0)
  //       node 2: (1,1)
  //       node 3: (0,1)
  //
  //       Nodes for 2nd order quads are as follows
  //       node 0: (0,0)
  //       node 1: (1,0)
  //       node 2: (1,1)
  //       node 3: (0,1)
  //       node 4: (1/2,0)
  //       node 5: (1,1/2)
  //       node 6: (1/2,1)
  //       node 7: (0,1/2)
  //       node 8: (1/2,1/2)
  if (p == 1)
  {
    // TODO: implement derivative of 1st order shapes for quads wrt ip.x
    dshape(0,0) = 0.;
    dshape(1,0) = 0.;
    dshape(2,0) = 0.;
    dshape(3,0) = 0.;

    // TODO: implement derivative of 1st order shapes for quads wrt ip.y
    dshape(0,1) = 0.;
    dshape(1,1) = 0.;
    dshape(2,1) = 0.;
    dshape(3,1) = 0.;
  }
  else if (p == 2)
  {
    // TODO: implement derivative of 2nd order shapes for quads wrt ip.x
    dshape(0,0) = 0.;
    dshape(1,0) = 0.;
    dshape(2,0) = 0.;
    dshape(3,0) = 0.;
    dshape(4,0) = 0.;
    dshape(5,0) = 0.;
    dshape(6,0) = 0.;
    dshape(7,0) = 0.;
    dshape(8,0) = 0.;

    // TODO: implement derivative of 2nd order shapes for quads wrt ip.y
    dshape(0,1) = 0.;
    dshape(1,1) = 0.;
    dshape(2,1) = 0.;
    dshape(3,1) = 0.;
    dshape(4,1) = 0.;
    dshape(5,1) = 0.;
    dshape(6,1) = 0.;
    dshape(7,1) = 0.;
    dshape(8,1) = 0.;
  }
  else
  {
    MFEM_VERIFY(0, "unimplemented order");
  }

}


// LG_TriangleElement implementation
LG_TriangleElement::LG_TriangleElement(
    const int p,
    const int btype)
   : NodalFiniteElement(2, Geometry::TRIANGLE, ((p + 1)*(p + 2))/2, p,
                        FunctionSpace::Pk)
{
  const double *cp = poly1d.ClosedPoints(p, VerifyNodal(VerifyClosed(btype)));

  // vertices
  Nodes.IntPoint(0).Set2(cp[0], cp[0]);
  Nodes.IntPoint(1).Set2(cp[p], cp[0]);
  Nodes.IntPoint(2).Set2(cp[0], cp[p]);

  // edges
  int o = 3;
  for (int i = 1; i < p; i++)
  {
    Nodes.IntPoint(o++).Set2(cp[i], cp[0]);
  }
  for (int i = 1; i < p; i++)
  {
    Nodes.IntPoint(o++).Set2(cp[p-i], cp[i]);
  }
  for (int i = 1; i < p; i++)
  {
    Nodes.IntPoint(o++).Set2(cp[0], cp[p-i]);
  }

  // interior
  for (int j = 1; j < p; j++)
    for (int i = 1; i + j < p; i++)
    {
      const double w = cp[i] + cp[j] + cp[p-i-j];
      Nodes.IntPoint(o++).Set2(cp[i]/w, cp[j]/w);
    }
}

void LG_TriangleElement::CalcShape(
    const IntegrationPoint &ip,
    Vector &shape) const
{
  // get the order from the class member "order"
  const int p = order;
  // Note: in all cases boundary nodes are the first two and then the interior nodes
  if (p == 1)
  {
    // TODO: implement 1st order shapes for triangles
    // Note the parametric coordinates are ip.x and ip.y
    shape(0) = 0.;
    shape(1) = 0.;
    shape(2) = 0.;
  }
  else if (p == 2)
  {
    // TODO: implement 2nd order shapes for triangles
    shape(0) = 0.;
    shape(1) = 0.;
    shape(2) = 0.;
    shape(3) = 0.;
    shape(4) = 0.;
    shape(5) = 0.;
  }
  else
  {
    MFEM_VERIFY(0, "unimplemented order");
  }
}

void LG_TriangleElement::CalcDShape(
    const IntegrationPoint &ip,
    DenseMatrix &dshape) const
{
  // get the order from the class member "order"
  const int p = order;
  // Note: in all cases boundary nodes are the first two and then the interior nodes
  if (p == 1)
  {
    // TODO: implement first derivative of 1st order shapes with respect to ip.x for triangles
    dshape(0,0) = 0.;
    dshape(1,0) = 0.;
    dshape(2,0) = 0.;

    // TODO: implement first derivative of 1st order shapes with respect to ip.y for triangles
    dshape(0,1) = 0.;
    dshape(1,1) = 0.;
    dshape(2,1) = 0.;
  }
  else if (p == 2)
  {
    // TODO: implement first derivative of 2nd order shapes with respect to ip.x for triangles
    dshape(0,0) = 0.;
    dshape(1,0) = 0.;
    dshape(2,0) = 0.;
    dshape(3,0) = 0.;
    dshape(4,0) = 0.;
    dshape(5,0) = 0.;

    // TODO: implement first derivative of 2nd order shapes with respect to ip.y for triangles
    dshape(0,1) = 0.;
    dshape(1,1) = 0.;
    dshape(2,1) = 0.;
    dshape(3,1) = 0.;
    dshape(4,1) = 0.;
    dshape(5,1) = 0.;
  }
  else
  {
    MFEM_VERIFY(0, "unimplemented order");
  }
}


// LG_FECollection implementation
LG_FECollection::LG_FECollection(const int p, const int dim, const int btype)
  : FiniteElementCollection(p)
{
  MFEM_VERIFY(p >= 1, "LG_FECollection requires order >= 1.");
  MFEM_VERIFY(p <  3, "LG_FECollection requires order <  3.");
  MFEM_VERIFY(dim == 2, "LG_FECollection requires dim == 2");

  const int pm1 = p - 1, pm2 = pm1 - 1, pm3 = pm2 - 1, pm4 = pm3 - 1;

  int pt_type = BasisType::GetQuadrature1D(btype);
  b_type = BasisType::Check(btype);
  MFEM_VERIFY(btype == BasisType::GaussLobatto,
      "LG_FECollection requires btype == BasisType::GaussLobatto");

  snprintf(lg_name, 32, "LG_%dD_P%d", dim, p);

  for (int g = 0; g < Geometry::NumGeom; g++)
  {
    LG_dof[g] = 0;
    LG_Elements[g] = NULL;
  }
  for (int i = 0; i < 2; i++)
  {
    SegDofOrd[i] = NULL;
  }
  for (int i = 0; i < 6; i++)
  {
    TriDofOrd[i] = NULL;
  }
  for (int i = 0; i < 8; i++)
  {
    QuadDofOrd[i] = NULL;
  }
  for (int i = 0; i < 24; i++)
  {
    TetDofOrd[i] = NULL;
  }

  LG_dof[Geometry::POINT] = 1;
  LG_Elements[Geometry::POINT] = new PointFiniteElement;

  if (dim >= 1)
  {
    LG_dof[Geometry::SEGMENT] = pm1;
    LG_Elements[Geometry::SEGMENT] = new LG_SegmentElement(p, btype);

    SegDofOrd[0] = new int[2*pm1];
    SegDofOrd[1] = SegDofOrd[0] + pm1;
    for (int i = 0; i < pm1; i++)
    {
      SegDofOrd[0][i] = i;
      SegDofOrd[1][i] = pm2 - i;
    }
  }

  if (dim >= 2)
  {
    LG_dof[Geometry::TRIANGLE] = (pm1*pm2)/2;
    LG_dof[Geometry::SQUARE] = pm1*pm1;
    LG_Elements[Geometry::TRIANGLE] = new LG_TriangleElement(p, btype);
    LG_Elements[Geometry::SQUARE] = new LG_QuadrilateralElement(p, btype);

    const int &TriDof = LG_dof[Geometry::TRIANGLE];
    TriDofOrd[0] = new int[6*TriDof];
    for (int i = 1; i < 6; i++)
    {
      TriDofOrd[i] = TriDofOrd[i-1] + TriDof;
    }
    // see Mesh::GetTriOrientation in mesh/mesh.cpp
    for (int j = 0; j < pm2; j++)
    {
      for (int i = 0; i + j < pm2; i++)
      {
        int o = TriDof - ((pm1 - j)*(pm2 - j))/2 + i;
        int k = pm3 - j - i;
        TriDofOrd[0][o] = o;  // (0,1,2)
        TriDofOrd[1][o] = TriDof - ((pm1-j)*(pm2-j))/2 + k;  // (1,0,2)
        TriDofOrd[2][o] = TriDof - ((pm1-i)*(pm2-i))/2 + k;  // (2,0,1)
        TriDofOrd[3][o] = TriDof - ((pm1-k)*(pm2-k))/2 + i;  // (2,1,0)
        TriDofOrd[4][o] = TriDof - ((pm1-k)*(pm2-k))/2 + j;  // (1,2,0)
        TriDofOrd[5][o] = TriDof - ((pm1-i)*(pm2-i))/2 + j;  // (0,2,1)
      }
    }

    const int &QuadDof = LG_dof[Geometry::SQUARE];
    QuadDofOrd[0] = new int[8*QuadDof];
    for (int i = 1; i < 8; i++)
    {
      QuadDofOrd[i] = QuadDofOrd[i-1] + QuadDof;
    }

    if (b_type == BasisType::Serendipity)
    {
      MFEM_VERIFY(0, "unimplemented quad type");
    }
    else // not serendipity
    {
      for (int j = 0; j < pm1; j++)
      {
        for (int i = 0; i < pm1; i++)
        {
          int o = i + j*pm1;
          QuadDofOrd[0][o] = i + j*pm1;  // (0,1,2,3)
          QuadDofOrd[1][o] = j + i*pm1;  // (0,3,2,1)
          QuadDofOrd[2][o] = j + (pm2 - i)*pm1;  // (1,2,3,0)
          QuadDofOrd[3][o] = (pm2 - i) + j*pm1;  // (1,0,3,2)
          QuadDofOrd[4][o] = (pm2 - i) + (pm2 - j)*pm1;  // (2,3,0,1)
          QuadDofOrd[5][o] = (pm2 - j) + (pm2 - i)*pm1;  // (2,1,0,3)
          QuadDofOrd[6][o] = (pm2 - j) + i*pm1;  // (3,0,1,2)
          QuadDofOrd[7][o] = i + (pm2 - j)*pm1;  // (3,2,1,0)
        }
      }
    }
  }
}

const int *LG_FECollection::DofOrderForOrientation(
    Geometry::Type GeomType,
    int Or) const
{
  if (GeomType == Geometry::SEGMENT)
  {
    return (Or > 0) ? SegDofOrd[0] : SegDofOrd[1];
  }
  else if (GeomType == Geometry::TRIANGLE)
  {
    return TriDofOrd[Or%6];
  }
  else if (GeomType == Geometry::SQUARE)
  {
    return QuadDofOrd[Or%8];
  }
  return NULL;
}

FiniteElementCollection *LG_FECollection::GetTraceCollection() const
{
  return NULL;
  /* int p = LG_dof[Geometry::SEGMENT] + 1; */
  /* int dim = -1; */
  /* if (!strncmp(lg_name, "LG_", 3)) */
  /* { */
  /*    dim = atoi(lg_name + 3); */
  /* } */
  /* /1* else if (!strncmp(h1_name, "H1Pos_", 6)) *1/ */
  /* /1* { *1/ */
  /* /1*    dim = atoi(h1_name + 6); *1/ */
  /* /1* } *1/ */
  /* /1* else if (!strncmp(h1_name, "H1@", 3)) *1/ */
  /* /1* { *1/ */
  /* /1*    dim = atoi(h1_name + 5); *1/ */
  /* /1* } *1/ */
  /* return (dim < 0) ? NULL : new LG_Trace_FECollection(p, dim, b_type); */
}

const int *LG_FECollection::GetDofMap(Geometry::Type GeomType) const
{
  const int *dof_map = NULL;
  const FiniteElement *fe = LG_Elements[GeomType];
  switch (GeomType)
  {
     case Geometry::SEGMENT:
     case Geometry::SQUARE:
     case Geometry::CUBE:
     default:
        MFEM_ABORT("Geometry type " << Geometry::Name[GeomType] << " is not "
                   "implemented");
        // The "Cartesian" ordering for other geometries is defined by the
        // class GeometryRefiner.
  }
  return dof_map;
}

LG_FECollection::~LG_FECollection()
{
   delete [] SegDofOrd[0];
   delete [] TriDofOrd[0];
   delete [] QuadDofOrd[0];
   delete [] TetDofOrd[0];
   for (int g = 0; g < Geometry::NumGeom; g++)
   {
      delete LG_Elements[g];
   }
}

#endif
