import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;

public class SmallPolygonsGrader {
    public static void main(String[] args) {
        try {
            String exec = args[0];
            Long seed = Long.valueOf(args[1]);
            Game game = new Game(seed);

            long limit = System.currentTimeMillis()+10000;
            Connector con = new Connector(exec);
            String[] ret = con.passIO(game);
            if (System.currentTimeMillis() > limit) {
                throw new RuntimeException("time up!");
            }
            double score = game.calcScore(ret);
            System.out.println(score);
        } catch (Exception e) {
            // something wrong. the score will be zero.
            System.out.println(0);
        }
    }

    public static class Connector {
        Runtime rt;
        Process proc;
        OutputStream os;
        BufferedReader br;

        public Connector(String exec) throws IOException {
            rt = Runtime.getRuntime();
            proc = rt.exec(exec);
            os = proc.getOutputStream();
            br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        }

        public String[] passIO(Game game) throws IOException {
            StringBuffer sb = new StringBuffer();
            sb.append(game.pointsPar.length).append('\n');
            for (int i = 0; i < game.pointsPar.length; ++i)
                sb.append(game.pointsPar[i]).append('\n');
            sb.append(game.N).append('\n');
            os.write(sb.toString().getBytes());
            os.flush();

            int nret = Integer.parseInt(br.readLine());
            String[] ret = new String[nret];
            for (int i = 0; i < nret; ++i) {
                ret[i] = br.readLine();
            }
            return ret;
        }
    }


    // ------------- class Point ------------------------------
    static class Pnt {
        public int x, y;

        public Pnt(int x1, int y1) {
            x = x1;
            y = y1;
        }

        public boolean equals(Pnt other) {
            return (x == other.x && y == other.y);
        }
    }

    // ------------- class G2D --------------------------------
    static class G2D {
        public static Pnt substr(Pnt p1, Pnt p2) {
            return new Pnt(p1.x-p2.x, p1.y-p2.y);
        }

        public static double norm(Pnt p) {
            return Math.sqrt(p.x * p.x+p.y * p.y);
        }

        public static int norm2(Pnt p) {
            return (p.x * p.x+p.y * p.y);
        }

        public static int dot(Pnt p1, Pnt p2) {
            return p1.x * p2.x+p1.y * p2.y;
        }

        public static int cross(Pnt p1, Pnt p2) {
            return p1.x * p2.y-p1.y * p2.x;
        }

        public static double dist(Pnt p1, Pnt p2) {
            return norm(substr(p1, p2));
        }

        public static int dist2(Pnt p1, Pnt p2) {
            return norm2(substr(p1, p2));
        }
    }

    // ------------- class Edge ------------------------------
    static class Edge {
        public Pnt p1, p2, vect;    //vector p1 -> p2
        public double norm;

        public Edge() {
        }

        ;

        public Edge(Pnt p1n, Pnt p2n) {
            p1 = p1n;
            p2 = p2n;
            vect = G2D.substr(p2, p1);
            norm = G2D.norm(vect);
        }

        public Edge(int x1, int y1, int x2, int y2) {
            p1 = new Pnt(x1, y1);
            p2 = new Pnt(x2, y2);
            vect = G2D.substr(p2, p1);
            norm = G2D.norm(vect);
        }

        boolean eq(double a, double b) {
            return Math.abs(a-b) < 1e-9;
        }

        // ---------------------------------------------------
        public boolean intersect(Edge other) {
            //do edges "this" and "other" intersect?
            if (Math.min(p1.x, p2.x) > Math.max(other.p1.x, other.p2.x)) return false;
            if (Math.max(p1.x, p2.x) < Math.min(other.p1.x, other.p2.x)) return false;
            if (Math.min(p1.y, p2.y) > Math.max(other.p1.y, other.p2.y)) return false;
            if (Math.max(p1.y, p2.y) < Math.min(other.p1.y, other.p2.y)) return false;

            int den = other.vect.y * vect.x-other.vect.x * vect.y;
            int num1 = other.vect.x * (p1.y-other.p1.y)-other.vect.y * (p1.x-other.p1.x);
            int num2 = vect.x * (p1.y-other.p1.y)-vect.y * (p1.x-other.p1.x);

            //parallel edges
            if (den == 0) {
                if (Math.min(other.dist2(this), dist2(other)) > 0)
                    return false;
                //on the same line - "not intersect" only if one of the vertices is common,
                //and the other doesn't belong to the line
                if ((this.p1 == other.p1 && eq(G2D.dist(this.p2, other.p2), this.norm+other.norm)) ||
                        (this.p1 == other.p2 && eq(G2D.dist(this.p2, other.p1), this.norm+other.norm)) ||
                        (this.p2 == other.p1 && eq(G2D.dist(this.p1, other.p2), this.norm+other.norm)) ||
                        (this.p2 == other.p2 && eq(G2D.dist(this.p1, other.p1), this.norm+other.norm)))
                    return false;
                return true;
            }

            //common vertices
            if (this.p1 == other.p1 || this.p1 == other.p2 || this.p2 == other.p1 || this.p2 == other.p2)
                return false;

            double u1 = num1 * 1. / den;
            double u2 = num2 * 1. / den;
            if (u1 < 0 || u1 > 1 || u2 < 0 || u2 > 1)
                return false;
            return true;
        }

        // ---------------------------------------------------
        public double dist(Pnt p) {
            //distance from p to the edge
            if (G2D.dot(vect, G2D.substr(p, p1)) <= 0)
                return G2D.dist(p, p1);         //from p to p1
            if (G2D.dot(vect, G2D.substr(p, p2)) >= 0)
                return G2D.dist(p, p2);         //from p to p2
            //distance to the line itself
            return Math.abs(-vect.y * p.x+vect.x * p.y+p1.x * p2.y-p1.y * p2.x) / norm;
        }

        // ---------------------------------------------------
        public double dist2(Edge other) {
            //distance from the closest of the endpoints of "other" to "this"
            return Math.min(dist(other.p1), dist(other.p2));
        }
    }

    public static class Game {
        final int SZ = 700;             // field size
        long seed;
        int NP, N, Npoly;               // number of points given, max number of polygons and number of polygons selected
        Pnt[] p;                        // coordinates of points (fixed)
        int[] pointsPar;                // coordinates of points (as an array parameter)
        int[][] polys;                  // indices of points which form polygons
        int[] polysVert;                // number of vertices in each poly
        boolean valid[];
        int[] used;                     // which poly uses this point?

        public Game(long seed) {
            this.seed = seed;
            SecureRandom rnd;
            try {
                rnd = SecureRandom.getInstance("SHA1PRNG");
            } catch (NoSuchAlgorithmException e) {
                throw new RuntimeException(e);
            }
            rnd.setSeed(seed);

            if (seed == 1) {
                NP = 10;
            } else {
                int testSize = rnd.nextInt(3);
                if (testSize == 0) {
                    NP = rnd.nextInt(80)+20;
                } else if (testSize == 1) {
                    NP = rnd.nextInt(400)+100;
                } else {
                    NP = rnd.nextInt(1001)+500;
                }
            }

            p = new Pnt[NP];
            used = new int[NP];
            Arrays.fill(used, -2);

            // generate the points
            boolean ok;
            for (int i = 0; i < NP; ++i) {
                do {
                    p[i] = new Pnt(rnd.nextInt(SZ), rnd.nextInt(SZ));
                    ok = true;
                    for (int j = 0; j < i && ok; ++j) {
                        if (p[i].equals(p[j])) {
                            ok = false;
                        }
                    }
                }
                while (!ok);
            }

            // convert points to parameter array
            pointsPar = new int[2 * NP];
            for (int i = 0; i < NP; ++i) {
                pointsPar[2 * i] = p[i].x;
                pointsPar[2 * i+1] = p[i].y;
            }
            N = rnd.nextInt(19)+2;
            if (seed == 1) {
                N = 3;
            }
        }

        public double calcScore(String[] ret) {
            Npoly = ret.length;
            int n = Npoly;
            polys = new int[n][];
            polysVert = new int[n];
            valid = new boolean[n];
            for (int i = 0; i < Npoly; ++i) {
                // parse the string into the polygon
                try {
                    String[] st = ret[i].split(" ");
                    int nv = st.length;
                    polys[i] = new int[nv];
                    polysVert[i] = nv;
                    for (int j = 0; j < nv; ++j) {
                        polys[i][j] = Integer.parseInt(st[j]);
                        // check whether this point already was used
                        if (used[polys[i][j]] > -2) {
                            return 0;
                        } else {
                            used[polys[i][j]] = i;
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    return 0;
                }
                if (!validatePoly(polys[i], polysVert[i])) {
                    return 0;
                }
                valid[i] = true;
            }

            if (Npoly > N) {
                return 0;
            }

            // 2. each point is used by one of polygons (no polygons using same point checked earlier)
            for (int i = 0; i < used.length; ++i) {
                if (used[i] == -2) {
                    return 0;
                }
            }

            // 3. each polygon is valid on its own
            for (int i = 0; i < polys.length; ++i) {
                if (!valid[i]) {
                    return 0;
                }
            }

            // 4. no two polygons intersect
            // for an intersection to be detected, polygons have to have at least 2 vertices, so deleted polygons in manual mode have no effect
            for (int i = 0; i < polys.length; ++i) {
                for (int j = 0; j < polysVert[i]; ++j) {
                    for (int k = i+1; k < polys.length; ++k) {
                        for (int l = 0; l < polysVert[k]; ++l) {
                            // check intersection of edge j..j+1 of polygon i and edge l..l+1 of polygon k
                            Edge e1 = new Edge(p[polys[i][j]], p[polys[i][(j+1) % polysVert[i]]]);
                            Edge e2 = new Edge(p[polys[k][l]], p[polys[k][(l+1) % polysVert[k]]]);
                            if (e1.intersect(e2)) {
                                return 0;
                            }
                        }
                    }
                }
            }

            // now, if all are valid, score is always non-0
            double score = 0;
            for (int i = 0; i < polys.length; ++i) {
                score += area(polys[i], polysVert[i]);
            }
            return score;
        }

        boolean validatePoly(int[] poly, int n) {
            // check that the polygon satisfies all individual conditions
            if (n < 3) {
                return false;
            }

            // simple polygon: no self-intersections except in common vertices of edges
            for (int i = 0; i < n; ++i)
                for (int j = i+1; j < n; ++j) {
                    // check intersection of i..i+1 and j..j+1
                    Edge e1 = new Edge(p[poly[i]], p[poly[(i+1) % n]]);
                    Edge e2 = new Edge(p[poly[j]], p[poly[(j+1) % n]]);
                    if (e1.intersect(e2)) {
                        return false;
                    }
                }
            return true;
        }

        // ---------------------------------------------------
        double area(int[] poly, int n) {
            // trapezium method
            double s = 0;
            for (int i = 0; i < n; i++)
                s += (p[poly[(i+1) % n]].y+p[poly[i]].y) * (p[poly[(i+1) % n]].x-p[poly[i]].x) / 2.0;
            return Math.abs(s);
        }
    }
}
