from siman.geo import xred2xcart


def _endpoint_cart(st):
    return st.xcart if st.xcart else xred2xcart(st.xred, st.rprimd)
