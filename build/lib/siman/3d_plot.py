if 0:
    from impurity import find_pores
    from scipy.spatial import ConvexHull

    st1 = calc['V249vTiV','8e4ns',1].end
    st2 = calc['V249vCrV','8e4ns',1].end
    # x1 = [ 4.31371359,  5.77865732,  5.71961224]
    x2 = [ 4.39062872,  6.04561325,  5.98035926]
    # st_pores = find_pores(st1, r_matrix = 0.5, r_impurity = 1.5, fine = 2, calctype = 'all_pores')
    # print (st_pores.xcart)

    sur = local_surrounding(x2, st2, n_neighbours = 13, control = 'atoms', periodic  = True)
    # write_xyz(st2, analysis = 'imp_surrounding' , show_around_x = x2, nnumber = 13)
    print(sur[2])
    # sys.exit()

    vacancies = []
    for st in st1, st2:
        points = [st.xcart[i] for i in sur[2]]


        points = np.array(points)
        center = sum(points)/len(points)
        # print (center)
        points = points - center

        hull = ConvexHull(points)
        polygons = []
        for simpl in hull.simplices:
            polygons.append([points[i] for i in simpl])
        # print (polygons)
        vacancies.append(polygons)
        print(hull.volume)
        # write_xyz(st, analysis = 'imp_surrounding' , show_around_x = x_vac, nnumber = 14)
    # write_xyz(st,)
    vac1 = vacancies[0]
    vac2 = vacancies[1]+np.array([0,6,0])
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import mpl_toolkits.mplot3d as a3
    import matplotlib.pyplot as plt
    import pylab as pl
    import scipy as sp

    import numpy
    from mpl_toolkits.mplot3d import proj3d
    def orthogonal_proj(zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        return numpy.array([[1,0,0,0],
                            [0,1,0,0],
                            [0,0,a,b],
                            [0,0,-0.0001,zback]])
    proj3d.persp_transformation = orthogonal_proj



    ax = a3.Axes3D(pl.figure())
    # fig = plt.figure()
    # ax = Axes3D(fig)
    # tri = ax.add_collection3d(Poly3DCollection(vac1) )
    if 1:
        tri = a3.art3d.Poly3DCollection(vac1)
        tri.set_color('r')
        tri.set_edgecolor('k')
        tri.set_alpha(0.5)
        ax.add_collection3d(tri)
    if 1:
        tri = a3.art3d.Poly3DCollection(vac2)
        tri.set_color('b')
        tri.set_edgecolor('k')
        tri.set_alpha(0.5)
        ax.add_collection3d(tri)
    # tri = ax.add_collection3d(Poly3DCollection(vac2) )
    pl.show()
    # pl.savefig('poly.png')
