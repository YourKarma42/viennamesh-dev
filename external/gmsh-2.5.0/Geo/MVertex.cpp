// Gmsh - Copyright (C) 1997-2010 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#include <string.h>
#include <math.h>
#include "MVertex.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
#include "GFaceCompound.h"
#include "GmshMessage.h"
#include "StringUtils.h"

int MVertex::_globalNum = 0;
double MVertexLessThanLexicographic::tolerance = 1.e-6;

bool MVertexLessThanLexicographic::operator()(const MVertex *v1, const MVertex *v2) const
{
  if(v1->x() - v2->x() >  tolerance) return true;
  if(v1->x() - v2->x() < -tolerance) return false;
  if(v1->y() - v2->y() >  tolerance) return true;
  if(v1->y() - v2->y() < -tolerance) return false;
  if(v1->z() - v2->z() >  tolerance) return true;
  return false;
}

double angle3Vertices(MVertex *p1, MVertex *p2, MVertex *p3)
{
  SVector3 a(p1->x() - p2->x(), p1->y() - p2->y(), p1->z() - p2->z());
  SVector3 b(p3->x() - p2->x(), p3->y() - p2->y(), p3->z() - p2->z());
  SVector3 c = crossprod(a, b);
  double sinA = c.norm();
  double cosA = dot(a, b);
  return atan2 (sinA, cosA);  
}

MVertex::MVertex(double x, double y, double z, GEntity *ge, int num)
  : _visible(1), _order(1), _x(x), _y(y), _z(z), _ge(ge)
{
#pragma omp critical
  {
    if(num){
      _num = num;
      _globalNum = std::max(_globalNum, _num);
    }
    else{
      _num = ++_globalNum;
    }
    _index = num;
  }
}

void MVertex::forceNum(int num)
{ 
#pragma omp critical
  {
    _num = num; 
    _globalNum = std::max(_globalNum, _num);
  }
}

void MVertex::writeMSH(FILE *fp, bool binary, bool saveParametric, double scalingFactor)
{
  if(_index < 0) return; // negative index vertices are never saved

  int myDim = 0, myTag = 0;
  if(saveParametric){
    if(onWhat()){
      myDim = onWhat()->dim(); 
      myTag = onWhat()->tag();
    }
    else
      saveParametric = false;
  }

  if(!binary){
    if(!saveParametric)
      fprintf(fp, "%d %.16g %.16g %.16g\n", _index, x() * scalingFactor, 
              y() * scalingFactor, z() * scalingFactor);      
    else
      fprintf(fp, "%d %.16g %.16g %.16g %d %d", _index, x() * scalingFactor, 
              y() * scalingFactor, z() * scalingFactor, myDim, myTag);
  }
  else{
    fwrite(&_index, sizeof(int), 1, fp);
    double data[3] = {x() * scalingFactor, y() * scalingFactor, z() * scalingFactor};
    fwrite(data, sizeof(double), 3, fp);
    if(saveParametric){
      fwrite(&myDim, sizeof(int), 1, fp);
      fwrite(&myTag, sizeof(int), 1, fp);
    }
  }

  if(saveParametric){
    if(myDim == 1){
      double _u;
      getParameter(0, _u);
      if(!binary)
        fprintf(fp, " %.16g\n", _u);        
      else
        fwrite(&_u, sizeof(double), 1, fp);
    }
    else if (myDim == 2){
      double _u, _v;
      getParameter(0, _u);
      getParameter(1, _v);
      if(!binary)
        fprintf(fp, " %.16g %.16g\n", _u, _v);
      else{
        fwrite(&_u, sizeof(double), 1, fp);
        fwrite(&_v, sizeof(double), 1, fp);
      }
    }
    else
      if(!binary)
        fprintf(fp, "\n");          
  }
}

void MVertex::writePLY2(FILE *fp)
{
  if(_index < 0) return; // negative index vertices are never saved

  fprintf(fp, "%.16g %.16g %.16g\n", x(), y(), z());
}

void MVertex::writeVRML(FILE *fp, double scalingFactor)
{
  if(_index < 0) return; // negative index vertices are never saved

  fprintf(fp, "%.16g %.16g %.16g,\n",
          x() * scalingFactor, y() * scalingFactor, z() * scalingFactor);
}

void MVertex::writeUNV(FILE *fp, double scalingFactor)
{
  if(_index < 0) return; // negative index vertices are never saved

  int coord_sys = 1;
  int displacement_coord_sys = 1;
  int color = 11;
  fprintf(fp, "%10d%10d%10d%10d\n", _index, coord_sys, displacement_coord_sys, color);
  // hack to print the numbers with "D+XX" exponents
  char tmp[128];
  sprintf(tmp, "%25.16E%25.16E%25.16E\n", x() * scalingFactor, 
          y() * scalingFactor, z() * scalingFactor);
  for(unsigned int i = 0; i < strlen(tmp); i++) if(tmp[i] == 'E') tmp[i] = 'D';
  fprintf(fp, "%s", tmp);
}

void MVertex::writeVTK(FILE *fp, bool binary, double scalingFactor, bool bigEndian)
{
  if(_index < 0) return; // negative index vertices are never saved

  if(binary){
    double data[3] = {x() * scalingFactor, y() * scalingFactor, z() * scalingFactor};
    // VTK always expects big endian binary data
    if(!bigEndian) SwapBytes((char*)data, sizeof(double), 3);
    fwrite(data, sizeof(double), 3, fp);
  }
  else{
    fprintf(fp, "%.16g %.16g %.16g\n",
            x() * scalingFactor, y() * scalingFactor, z() * scalingFactor);
  }
}

void MVertex::writeMESH(FILE *fp, double scalingFactor)
{
  if(_index < 0) return; // negative index vertices are never saved

  fprintf(fp, " %20.14G      %20.14G      %20.14G      %d\n", 
          x() * scalingFactor, y() * scalingFactor, z() * scalingFactor, 
          _ge ? _ge->tag() : 0);
}

static void double_to_char8(double val, char *str)
{
  if(val >= 1.e6)
    sprintf(str, "%.2E", val);
  else if(val >= 1.e-3)
    sprintf(str, "%f", val);
  else if(val >= 0)
    sprintf(str, "%.2E", val);
  else if(val >= -1.e-3)
    sprintf(str, "%.1E", val);
  else if(val >= -1.e6)
    sprintf(str, "%f", val);
  else
    sprintf(str, "%.1E", val);

#if defined(WIN32)
  // Windows uses 3 digits in the exponent (which apparently does not
  // conform with ANSI C): remove the extra 0
  if(strlen(str) == 9 && (str[4] == 'E' || str[5] == 'E')){
    str[6] = str[7];
    str[7] = str[8];
  }
#endif

  str[8] = '\0';
}

void MVertex::writeBDF(FILE *fp, int format, double scalingFactor)
{
  if(_index < 0) return; // negative index vertices are never saved

  char xs[17], ys[17], zs[17];
  double x1 = x() * scalingFactor;
  double y1 = y() * scalingFactor;
  double z1 = z() * scalingFactor;
  if(format == 0){
    // free field format (max 8 char per field, comma separated, 10 per line)
    double_to_char8(x1, xs); double_to_char8(y1, ys); double_to_char8(z1, zs);
    fprintf(fp, "GRID,%d,%d,%s,%s,%s\n", _index, 0, xs, ys, zs);
  }
  else if(format == 1){ 
    // small field format (8 char par field, 10 per line)
    double_to_char8(x1, xs); double_to_char8(y1, ys); double_to_char8(z1, zs);
    fprintf(fp, "GRID    %-8d%-8d%-8s%-8s%-8s\n", _index, 0, xs, ys, zs);
  }
  else{ 
    // large field format (8 char first/last field, 16 char middle, 6 per line)
    fprintf(fp, "GRID*   %-16d%-16d%-16.9G%-16.9G*N%-6d\n", _index, 0, x1, y1, _index);
    fprintf(fp, "*N%-6d%-16.9G\n", _index, z1);
  }
}

void MVertex::writeDIFF(FILE *fp, bool binary, double scalingFactor)
{
  if(_index < 0) return; // negative index vertices are never saved

  fprintf(fp, " %d ( %25.16E , %25.16E , %25.16E )",
          getIndex(), x() * scalingFactor, y() * scalingFactor, z() * scalingFactor);
}

std::set<MVertex*, MVertexLessThanLexicographic>::iterator 
MVertex::linearSearch(std::set<MVertex*, MVertexLessThanLexicographic> &pos)
{
  for(std::set<MVertex*, MVertexLessThanLexicographic>::iterator it = pos.begin();
      it != pos.end(); ++it)
    if(distance(*it) < MVertexLessThanLexicographic::tolerance) return it;
  return pos.end();
}

static void getAllParameters(MVertex *v, GFace *gf, std::vector<SPoint2> &params)
{
  params.clear();

  if (gf->geomType() == GEntity::CompoundSurface &&
      v->onWhat()->dim() < 2){
    GFaceCompound *gfc = (GFaceCompound*) gf;
    params.push_back(gfc->getCoordinates(v));
    return;
  }

  if(v->onWhat()->dim() == 0){
    GVertex *gv = (GVertex*)v->onWhat();
    std::list<GEdge*> ed = gv->edges();
    bool seam = false;
    for(std::list<GEdge*>::iterator it = ed.begin(); it != ed.end(); it++){
      if((*it)->isSeam(gf)) {
        Range<double> range = (*it)->parBounds(0);
        if (gv == (*it)->getBeginVertex()){
          params.push_back((*it)->reparamOnFace(gf, range.low(),-1));
          params.push_back((*it)->reparamOnFace(gf, range.low(), 1));
        }
        else if (gv == (*it)->getEndVertex()){
          params.push_back((*it)->reparamOnFace(gf, range.high(),-1));
          params.push_back((*it)->reparamOnFace(gf, range.high(), 1));
        }
        else{
          Msg::Warning("Strange!");
        }
        seam = true;
      }
    }
    if (!seam)
      params.push_back(gv->reparamOnFace(gf, 1));
  }
  else if(v->onWhat()->dim() == 1){
    GEdge *ge = (GEdge*)v->onWhat();
    double UU;
    v->getParameter(0, UU);
    if (UU == 0.0)
      UU = ge->parFromPoint(v->point());
    params.push_back(ge->reparamOnFace(gf, UU, 1));
    if(ge->isSeam(gf))
      params.push_back(ge->reparamOnFace(gf, UU, -1));
  }
  else{
    double UU, VV;
    if(v->onWhat() == gf && v->getParameter(0, UU) && v->getParameter(1, VV))
      params.push_back(SPoint2(UU, VV));
  }
}

bool reparamMeshEdgeOnFace(MVertex *v1, MVertex *v2, GFace *gf, 
                           SPoint2 &param1, SPoint2 &param2)
{
  std::vector<SPoint2> p1, p2;
  getAllParameters(v1, gf, p1);
  getAllParameters(v2, gf, p2);
  if (p1.size() == 1 && p2.size() == 1){
    param1 = p1[0];
    param2 = p2[0];
    return true;
  }
  else if (p1.size() == 1 && p2.size() == 2){
    double d1 = 
      (p1[0].x() - p2[0].x()) * (p1[0].x() - p2[0].x()) +
      (p1[0].x() - p2[0].y()) * (p1[0].y() - p2[0].y());
    double d2 = 
      (p1[0].x() - p2[1].x()) * (p1[0].x() - p2[1].x()) +
      (p1[0].x() - p2[1].y()) * (p1[0].y() - p2[1].y());
    param1 = p1[0];
    param2 = d2 < d1 ? p2[1] : p2[0];
    return true;
  }  
  else if (p2.size() == 1 && p1.size() == 2){
    double d1 = 
      (p2[0].x() - p1[0].x()) * (p2[0].x() - p1[0].x()) +
      (p2[0].x() - p1[0].y()) * (p2[0].y() - p1[0].y());
    double d2 = 
      (p2[0].x() - p1[1].x()) * (p2[0].x() - p1[1].x()) +
      (p2[0].x() - p1[1].y()) * (p2[0].y() - p1[1].y());
    param1 = d2 < d1 ? p1[1] : p1[0];
    param2 = p2[0];
    return true;
  }  
  else if(p1.size() > 1 && p2.size() > 1){
    param1 = p1[0];
    param2 = p2[0];
    // shout, both vertices are on seams
    return false;
  }
  else{
    // brute force!
    param1 = gf->parFromPoint(SPoint3(v1->x(), v1->y(), v1->z()));
    param2 = gf->parFromPoint(SPoint3(v2->x(), v2->y(), v2->z()));
    return true;
  }
}

bool reparamMeshVertexOnFace(const MVertex *v, const GFace *gf, SPoint2 &param,
                             bool onSurface)
{
  if (gf->geomType() == GEntity::CompoundSurface &&
      v->onWhat()->dim() < 2){
    GFaceCompound *gfc = (GFaceCompound*) gf;
    param = gfc->getCoordinates(const_cast<MVertex*>(v));
    return true;
  }

  if(v->onWhat()->geomType() == GEntity::DiscreteCurve ||        
     v->onWhat()->geomType() == GEntity::BoundaryLayerCurve){    
    param = gf->parFromPoint(SPoint3(v->x(), v->y(), v->z()), onSurface);
    return true;
  }

  if(v->onWhat()->dim() == 0){
    GVertex *gv = (GVertex*)v->onWhat();
    // hack for bug in periodic curves
    if (gv->getNativeType() == GEntity::GmshModel && gf->geomType() == GEntity::Plane)
      param = gf->parFromPoint(SPoint3(v->x(), v->y(), v->z()), onSurface);
    else
      param = gv->reparamOnFace(gf, 1);
    // shout, we could be on a seam
    std::list<GEdge*> ed = gv->edges();
    for(std::list<GEdge*>::iterator it = ed.begin(); it != ed.end(); it++)
      if((*it)->isSeam(gf)) return false;
  }
  else if(v->onWhat()->dim() == 1){
    GEdge *ge = (GEdge*)v->onWhat();
    double t;
    v->getParameter(0, t);
    param = ge->reparamOnFace(gf, t, 1);

    // shout, we are on a seam
    if(ge->isSeam(gf))
      return false;
  }
  else{
    double uu, vv;
    if(v->onWhat() == gf && v->getParameter(0, uu) && v->getParameter(1, vv)){
      param = SPoint2(uu, vv);
    }
    else {
      // brute force!
      param = gf->parFromPoint(SPoint3(v->x(), v->y(), v->z()), onSurface);
    }
  }
  return true;
}

bool reparamMeshVertexOnEdge(const MVertex *v, const GEdge *ge, double &param)
{
  param = 1.e6;
  Range<double> bounds = ge->parBounds(0);
  bool ok = true;
  if(ge->getBeginVertex() && ge->getBeginVertex()->mesh_vertices[0] == v)
    param = bounds.low();
  else if(ge->getEndVertex() && ge->getEndVertex()->mesh_vertices[0] == v)
    param = bounds.high();
  else
    ok = v->getParameter(0, param);

  if(!ok || param == 1.e6)
    param = ge->parFromPoint(SPoint3(v->x(), v->y(), v->z()));
  
  if(param < 1.e6) return true;
  return false;
}

#include "Bindings.h"

void MVertex::registerBindings(binding *b)
{
  classBinding *cb = b->addClass<MVertex>("MVertex");
  cb->setDescription("A mesh vertex.");
  methodBinding *cm;
  cm = cb->addMethod("getNum",&MVertex::getNum);
  cm->setDescription("Return the immutable vertex number.");
  //the cast is epxlicitely given because there are 2 MVertex::x function
  cm = cb->addMethod("x", (double (MVertex::*)() const) &MVertex::x);
  cm->setDescription("Return the x-coordinate.");
  cm = cb->addMethod("y", (double (MVertex::*)() const) &MVertex::y);
  cm->setDescription("Return the y-coordinate.");
  cm = cb->addMethod("z", (double (MVertex::*)() const) &MVertex::z);
  cm->setDescription("Return the z-coordinate.");
  cm = cb->addMethod("setXYZ", &MVertex::setXYZ);
  cm->setDescription("set the coordinates");
  cm->setArgNames("x", "y", "z",NULL);
  cm = cb->setConstructor<MVertex,double,double,double>();
  cm->setArgNames("x", "y", "z", NULL);
  cm->setDescription("Create a new mesh vertex at (x,y,z).");
  cm = cb->addMethod("getNum", &MVertex::getNum);
  cm->setDescription("return the invariant vertex id");
  cm = cb->addMethod("getPolynomialOrder", &MVertex::getPolynomialOrder);
  cm->setDescription("return the polynomial order of vertex");
  cm = cb->addMethod("setPolynomialOrder", &MVertex::setPolynomialOrder);
  cm->setDescription("assign the polynomial order of vertex");
  cm->setArgNames("order",NULL);
}