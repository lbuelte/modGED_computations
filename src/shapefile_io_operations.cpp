#include "../include/shapefile_io_operations.h"

using namespace boost::geometry;


void writeToShapeFile(MultiPolygon_boost polys, std::string path) {

    // Create handle
    SHPHandle shapefile = SHPCreate(path.c_str(), SHPT_POLYGON);
    DBFHandle dbfile = DBFCreate(path.c_str());

    // Create field for poly ID
    int field_face_id = DBFAddField(dbfile, "ID", FTInteger, 5, 0);

    int polyID = 0;

    for (const auto &poly: polys) {
        int partCount = 1 + boost::geometry::num_interior_rings(poly);

        // Collect vertices
        std::vector<double> x;
        std::vector<double> y;

        // Add exterior ring vertices
        for (const Point_boost &p: boost::geometry::exterior_ring(poly)) {
            x.push_back(p.a);
            y.push_back(p.b);
        }

        // Add interior ring vertices (holes)
        /*for (const auto& hole : boost::geometry::interior_rings(poly)) {
            for (const Point_boost& p : hole) {
                x.push_back(p.a);
                y.push_back(p.b);
            }
        }*/

        // Create polygon object
        SHPObject *shape = SHPCreateSimpleObject(SHPT_POLYGON, x.size(), x.data(), y.data(), nullptr);

        // Write shape into file
        int shape_id = SHPWriteObject(shapefile, -1, shape);

        // Set field to face ID
        DBFWriteIntegerAttribute(dbfile, shape_id, field_face_id, polyID++);

        // Memory management
        SHPDestroyObject(shape);
    }

    DBFClose(dbfile);
    SHPClose(shapefile);
}

void writeToShapeFile(MultiPolygon_boost polys, std::vector<int> clusters, std::string path) {

    // Create handle
    SHPHandle shapefile = SHPCreate(path.c_str(), SHPT_POLYGON);
    DBFHandle dbfile = DBFCreate(path.c_str());

    // Create field for poly ID
    int field_face_id = DBFAddField(dbfile, "ID", FTInteger, 5, 0);
    int field_cluster_id = DBFAddField(dbfile, "ClusterID", FTInteger, 10, 0);

    int polyID = 0;

    for (const auto &poly: polys) {
        int partCount = 1 + boost::geometry::num_interior_rings(poly);

        // Collect vertices
        std::vector<double> x;
        std::vector<double> y;

        // Add exterior ring vertices
        for (const Point_boost &p: boost::geometry::exterior_ring(poly)) {
            x.push_back(p.a);
            y.push_back(p.b);
        }

        // Add interior ring vertices (holes)
        /*for (const auto& hole : boost::geometry::interior_rings(poly)) {
            for (const Point_boost& p : hole) {
                x.push_back(p.a);
                y.push_back(p.b);
            }
        }*/

        // Create polygon object
        SHPObject *shape = SHPCreateSimpleObject(SHPT_POLYGON, x.size(), x.data(), y.data(), nullptr);

        // Write shape into file
        int shape_id = SHPWriteObject(shapefile, -1, shape);

        // Set field to face ID
        DBFWriteIntegerAttribute(dbfile, shape_id, field_face_id, polyID);

        // Set field to face ID
        DBFWriteIntegerAttribute(dbfile, shape_id, field_cluster_id, clusters[polyID]);

        // Memory management
        SHPDestroyObject(shape);

        polyID++;
    }

    DBFClose(dbfile);
    SHPClose(shapefile);
}



/*void ReadShapeFile(const std::string &filename, MultiPolygon_boost& multiPolygon) {
    try {
        SHPHandle handle = SHPOpen(filename.c_str(), "rb");
        if (handle->nFileSize <= 0) {
            throw std::string("File " + filename + " not found.");
        }

        int nShapeType, nEntities;
        double adfMinBound[4], adfMaxBound[4];
        SHPGetInfo(handle, &nEntities, &nShapeType, adfMinBound, adfMaxBound);

        //cout << "file has " << nEntities << " entities and shape type "<< nShapeType << endl;
        for (int i = 0; i < nEntities; i++) {
            SHPObject *psShape = SHPReadObject(handle, i);

            //cout << "read object, shapetype: " << psShape->nSHPType << endl;

            //Read Polys (only those without holes)
            if (psShape->nSHPType == SHPT_POLYGON && psShape->nParts == 1) {
                Polygon_boost polygon;

                double *x = psShape->padfX;
                double *y = psShape->padfY;
                //read polygon from file
                for (int v = psShape->nVertices - 1; v > 0; v--) {
                    Point_boost point = Point_boost(x[v], y[v]);

                    //add point to polygon
                    boost::geometry::append(polygon, point);

                }

                multiPolygon.push_back(polygon);
            }
            SHPDestroyObject(psShape);
        }
        SHPClose(handle);

    }
    catch (const std::string &s) {
        throw s;
    }

}*/

bool operator==(const Point_boost &lhs, const Point_boost &rhs) {
    return lhs.a == rhs.a && lhs.b == rhs.b;
}

void ReadShapeFile(const std::string &filename, MultiPolygon_boost &multiPolygon) {
    try {
        SHPHandle handle = SHPOpen(filename.c_str(), "rb");
        if (handle->nFileSize <= 0) {
            throw std::string("File " + filename + " not found.");
        }

        int nShapeType, nEntities;
        double adfMinBound[4], adfMaxBound[4];
        SHPGetInfo(handle, &nEntities, &nShapeType, adfMinBound, adfMaxBound);

        for (int i = 0; i < nEntities; i++) {
            SHPObject *psShape = SHPReadObject(handle, i);

            if (psShape->nSHPType == SHPT_POLYGON) {
                // Create a polygon
                Polygon_boost polygon;
                std::list<Polygon_boost> holes;
                std::vector<Point_boost> ring;

                // Read exterior ring
                for (int v = 0; v < psShape->nVertices; v++) {
                    Point_boost point = Point_boost(psShape->padfX[v], psShape->padfY[v]);
                    //boost::geometry::append(polygon, point);
                    ring.push_back(point);
                }

                // Read interior rings (holes)
                for (int part = 1; part < psShape->nParts; part++) {
                    int start = psShape->panPartStart[part];
                    int end = (part == psShape->nParts - 1) ? psShape->nVertices : psShape->panPartStart[part + 1];
                    auto &exteriorRing = boost::geometry::exterior_ring(polygon);
                    Polygon_boost hole = Polygon_boost();
                    for (int v = end - 1; v >= start; v--) {
                        const Point_boost &point = Point_boost(psShape->padfX[v], psShape->padfY[v]);
                        //auto it = std::count(boost::begin(polygon.outer()), boost::end(polygon.outer()), point);
                        //std::count(outer.begin(), outer.end(), point);


                        // If the point is found, erase it
                        if (std::find(exteriorRing.begin(), exteriorRing.end(), point) != exteriorRing.end()) {
                            auto rem = std::remove(exterior_ring(polygon).begin(), exterior_ring(polygon).end(), point);
                            exteriorRing.erase(rem);
                        }
                        if (std::find(ring.begin(), ring.end(), point) != ring.end()) {
                            ring.erase(std::find(ring.begin(), ring.end(), point));
                        }

                        boost::geometry::append(hole, point);
                    }

                    holes.push_back(hole);

                    //boost::geometry::difference(polygon, hole, result);
                    //boost::geometry::interior_rings(polygon).emplace_back(bg::exterior_ring(hole));
                }

                for (auto pp: ring) {
                    boost::geometry::append(polygon, pp);
                }

                for (auto hole: holes) {
                    boost::geometry::interior_rings(polygon).emplace_back(bg::exterior_ring(hole));
                }


                multiPolygon.push_back(polygon);
            }

            SHPDestroyObject(psShape);
        }
        SHPClose(handle);
    }
    catch (const std::string &s) {
        throw s;
    }
}


void writeToShapeFile(std::vector<Graph_boost> mas, std::string path, double angle = 1) {

    //create handle
    SHPHandle shapefile = SHPCreate(path.c_str(), SHPT_ARC);
    DBFHandle dbfile = DBFCreate(path.c_str());

    //create field for poly ID
    int field_label = DBFAddField(dbfile, "Label", FTInteger, 5, 0);
    int field_weight = DBFAddField(dbfile, "weight", FTDouble, 6, 0);

    int f = 0;

    for (const auto &g: mas) {
        //iterate through edges and add every edge to Shapefile
        auto edges = boost::edges(g);
        for (auto e = edges.first; e != edges.second; e++) {
            auto ed = boost::edge(e->m_source, e->m_target, g).first;

            if (g[ed].vanishingAngle >= angle) {

                //collect vertices
                double *x = new double[2];
                double *y = new double[2];

                x[0] = g[boost::source(*e, g)].p.a;
                y[0] = g[boost::source(*e, g)].p.b;
                x[1] = g[boost::target(*e, g)].p.a;
                y[1] = g[boost::target(*e, g)].p.b;

                //std::cout << g[ed].vanishingAngle << "\n";
                //std::cout<<x[1]<<"   "<<y[1]<<"\n";




                //create polygon object
                SHPObject *shape = SHPCreateSimpleObject(SHPT_ARC, 2, x, y, nullptr);

                //write shape into file
                int shape_id = SHPWriteObject(shapefile, -1, shape);

                //set field to face ID
                DBFWriteIntegerAttribute(dbfile, shape_id, field_label, g[*e].label);
                DBFWriteIntegerAttribute(dbfile, shape_id, field_weight, g[*e].weight);

                //memory management
                SHPDestroyObject(shape);
                delete[] x;
                delete[] y;
            }

        }
    }

    DBFClose(dbfile);
    SHPClose(shapefile);
}

