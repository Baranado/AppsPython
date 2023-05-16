import math
ety = 5;    ecu = -3;   phic = 0.65;    phit = 0.9; fatiga = 0.85


def FLEXCOMPBIAXGR(fc,  fy,  Es,     XYG,  XYS,  alfa,  c,     fcp=0.28):
    alfa = alfa * math.atan(1) / 45
    xygo = RangeToArray(XYG)
    xyso = RangeToArray(XYS)
    return FlexCompressGral(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def FLEXCOMPBIAXRN(fc,  fy,  Es,  b,      h,  nx,  ny,  dp,  ds,      de,  rec,  alfa,  c,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrR(b, h, ang)
    xyso = RebarRN(b, h, nx, ny, dp, ds, de, rec, ang)
    return FlexCompressGral(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def FLEXCOMPBIAXRS(fc,  fy,  Es,  b,      h,  dp,  ds,      de,  rec,  So,  alfa,  c,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrR(b, h, ang)
    xyso = RebarRS(b, h, dp, ds, de, rec, So, ang)
    return FlexCompressGral(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def flexcompbiaxrsu(fc,  fy,  Es,  b,      h,  dp,  ds,      de,  rec,  So,  alfa,  c,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrR(b, h, ang)
    xyso = RebarRS(b, h, dp, ds, de, rec, So, ang)
    return FlexCompressGralsu(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def FLEXCOMPBIAXT2(fc,  fy,  Es,  b,      h,  bw,  T,  b1,  t1,      dp,  ds,  de,  rec,  So,      alfa,  c,   fcp=0.28,       ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT2(b, h, bw, T, b1, t1, ang)
    xyso = RebarT2(b, h, bw, T, b1, t1, dp, ds, de, rec, So, ang)
    return FlexCompressGral(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def FLEXCOMPBIAXT1(fc,  fy,  Es,  b,      h,  bw,  T,  b1,      dp,  ds,  de,  rec,  So,      alfa,  c,   fcp=0.28,       ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT1(b, h, bw, T, b1, ang)
    xyso = RebarT1(b, h, bw, T, b1, dp, ds, de, rec, So, ang)
    return FlexCompressGral(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def FLEXCOMPBIAXT(fc,  fy,  Es,  b,      h,  bw,  T,      dp,  ds,  de,  rec,  So,      alfa,  c,   fcp=0.28,       ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT(b, h, bw, T, ang)
    xyso = RebarT(b, h, bw, T, dp, ds, de, rec, So, ang)
    return FlexCompressGral(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def FLEXCOMPBIAXL1(fc,  fy,  Es,  b,      h,  bw,  T,      dp,  ds,  de,  rec,  So,     alfa,  c,   fcp=0.28,      ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrL1(b, h, bw, T, ang)
    xyso = RebarL1(b, h, bw, T, dp, ds, de, rec, So, ang)
    return FlexCompressGral(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def FLEXCOMPBIAXL(fc,  fy,  Es,  b,      h,  T,      dp,  ds,  de,  rec,  So,      alfa,  c,   fcp=0.28,       ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrL(b, h, T, ang)
    xyso = RebarL(b, h, T, dp, ds, de, rec, So, ang)
    return FlexCompressGral(xygo, xyso, alfa, c, fc, fy, Es, fcp)


def XCOMPRESSGR(fc,  fy,  Es,  XYG,      XYS,  alfa,  P,   fcp=0.28):
    alfa = alfa * math.atan(1) / 45
    xygo = RangeToArray(XYG)
    xyso = RangeToArray(XYS)
    return xCompressGral(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def XCOMPRESSRN(fc,  fy,  Es,  b,      h,  nx,  ny,  dp,  ds,      de,  rec,  alfa,  P,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrR(b, h, ang)
    xyso = RebarRN(b, h, nx, ny, dp, ds, de, rec, ang)
    return xCompressGral(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def XCOMPRESSRS(fc,  fy,  Es,  b,      h,  dp,  ds,      de,  rec,  So,  alfa,  P,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrR(b, h, ang)
    xyso = RebarRS(b, h, dp, ds, de, rec, So, ang)
    return xCompressGralsu(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def xcompressrsu(fc,  fy,  Es,  b,      h,  dp,  ds,      de,  rec,  So,  alfa,  P,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrR(b, h, ang)
    xyso = RebarRS(b, h, dp, ds, de, rec, So, ang)
    return xCompressGralsu(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def XCOMPRESST2(fc,  fy,  Es,  b,      h,  bw,  T,  b1,  t1,  dp,  ds,      de,  rec,  So,  alfa,  P,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT2(b, h, bw, T, b1, t1, ang)
    xyso = RebarT2(b, h, bw, T, b1, t1, dp, ds, de, rec, So, ang)
    return xCompressGral(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def XCOMPRESST1(fc,  fy,  Es,  b,      h,  bw,  T,  b1,  dp,  ds,      de,  rec,  So,  alfa,  P,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT1(b, h, bw, T, b1, ang)
    xyso = RebarT1(b, h, bw, T, b1, dp, ds, de, rec, So, ang)
    return xCompressGral(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def XCOMPRESST(fc,  fy,  Es,  b,      h,  bw,  T,  dp,  ds,      de,  rec,  So,  alfa,  P,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT(b, h, bw, T, ang)
    xyso = RebarT(b, h, bw, T, dp, ds, de, rec, So, ang)
    return xCompressGral(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def XCOMPRESSL1(fc,  fy,  Es,  b,      h,  bw,  T,  dp,  ds,      de,  rec,  So,  alfa,  P,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrL1(b, h, bw, T, ang)
    xyso = RebarL1(b, h, bw, T, dp, ds, de, rec, So, ang)
    return xCompressGral(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def XCOMPRESSL(fc,  fy,  Es,  b,      h,  T,  dp,  ds,      de,  rec,  So,  alfa,  P,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrL(b, h, T, ang)
    xyso = RebarL(b, h, T, dp, ds, de, rec, So, ang)
    return xCompressGral(xygo(), xyso(), alfa, fc, fy, Es, P, fcp)


def DIAGRAFLEXCOMPGR(fc,  fy,  Es,      XYG,  XYS,  alfa,   fcp=0.28):
    alfa = alfa * math.atan(1) / 45
    xygo = RangeToArray(XYG)
    xyso = RangeToArray(XYS)
    return diagraFlexCompress(xygo, xyso, alfa, fc, fy, Es, fcp)


def DIAGRAFLEXCOMPRN(fc,  fy,  Es,  b,      h,  nx,  ny,  dp,  ds,      de,  rec,  alfa,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrR(b, h, ang)
    xyso = RebarRN(b, h, nx, ny, dp, ds, de, rec, ang)
    return diagraFlexCompress(xygo, xyso, alfa, fc, fy, Es, fcp)


def DIAGRAFLEXCOMPRS(fc,  fy,  Es,  b,      h,  dp,  ds,      de,  rec,  So,  alfa,       fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrR(b, h, ang)
    xyso = RebarRS(b, h, dp, ds, de, rec, So, ang)
    return diagraFlexCompress(xygo, xyso, alfa, fc, fy, Es, fcp)


def DIAGRAFLEXCOMPT2(fc,  fy,  Es,      b,  h,  bw,  T,  b1,      t1,  dp,  ds,  de,  rec,      So,  alfa,   fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT2(b, h, bw, T, b1, t1, ang)
    xyso = RebarT2(b, h, bw, T, b1, t1, dp, ds, de, rec, So, ang)
    return diagraFlexCompress(xygo, xyso, alfa, fc, fy, Es, fcp)


def DIAGRAFLEXCOMPT1(fc,  fy,  Es,      b,  h,  bw,  T,  b1,      dp,  ds,  de,  rec,      So,  alfa,   fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT1(b, h, bw, T, b1, ang)
    xyso = RebarT1(b, h, bw, T, b1, dp, ds, de, rec, So, ang)
    return diagraFlexCompress(xygo, xyso, alfa, fc, fy, Es, fcp)


def DIAGRAFLEXCOMPT(fc,  fy,  Es,      b,  h,  bw,  T,      dp,  ds,  de,  rec,      So,  alfa,   fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrT(b, h, bw, T, ang)
    xyso = RebarT(b, h, bw, T, dp, ds, de, rec, So, ang)
    return diagraFlexCompress(xygo, xyso, alfa, fc, fy, Es, fcp)


def DIAGRAFLEXCOMPL1(fc,  fy,  Es,      b,  h,  bw,  T,      dp,  ds,  de,  rec,      So,  alfa,   fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrL1(b, h, bw, T, ang)
    xyso = RebarL1(b, h, bw, T, dp, ds, de, rec, So, ang)
    return diagraFlexCompress(xygo, xyso, alfa, fc, fy, Es, fcp)


def DIAGRAFLEXCOMPL(fc,  fy,  Es,      b,  h,  T,     dp,  ds,  de,  rec,      So,  alfa,   fcp=0.28,   ang=0):
    alfa = alfa * math.atan(1) / 45
    ang = -ang * math.atan(1) / 45
    xygo = concrL(b, h, T, ang)
    xyso = RebarL(b, h, T, dp, ds, de, rec, So, ang)
    return diagraFlexCompress(xygo, xyso, alfa, fc, fy, Es, fcp)


def CONCREC(b,  h,   ang=0):
    ang = -ang * math.atan(1) / 45
    return concrR(b, h, ang)


def CONCTE2(b,  h,  bw,      T,  b1,  t1,       ang=0):
    ang = -ang * math.atan(1) / 45
    return concrT2(b, h, bw, T, b1, t1, ang)


def CONCTE1(b,  h,  bw,      T,  b1,   ang=0):
    ang = -ang * math.atan(1) / 45
    return concrT1(b, h, bw, T, b1, ang)


def CONCTEE(b,  h,  bw,      T,   ang=0):
    ang = -ang * math.atan(1) / 45
    return concrT(b, h, bw, T, ang)


def CONCEL1(b,  h,  bw,      T,   ang=0):
    ang = -ang * math.atan(1) / 45
    return concrL1(b, h, bw, T, ang)


def CONCELE(b,  h,      T,   ang=0):
    ang = -ang * math.atan(1) / 45
    return concrL(b, h, T, ang)


def REBARECN(b,  h,  nx,  ny,      dp,  ds,  de,      rec,   ang=0):
    ang = -ang * math.atan(1) / 45
    return RebarRN(b, h, nx, ny, dp, ds, de, rec, ang)


def REBARECS(b,  h,      dp,  ds,  de,      rec,  So,   ang=0):
    ang = -ang * math.atan(1) / 45
    return RebarRS(b, h, dp, ds, de, rec, So, ang)


def REBARTE2(b,  h,  bw,      T,  b1,  t1,      dp,  ds,  de,      rec,  So,   ang=0):
    ang = -ang * math.atan(1) / 45
    return RebarT2(b, h, bw, T, b1, t1, dp, ds, de, rec, So, ang)


def REBARTE1(b,  h,  bw,      T,  b1,  dp,  ds,      de,  rec,  So,   ang=0):
    ang = -ang * math.atan(1) / 45
    return RebarT1(b, h, bw, T, b1, dp, ds, de, rec, So, ang)


def REBARTEE(b,  h,  bw,      T,  dp,  ds,      de,  rec,  So,   ang=0):
    ang = -ang * math.atan(1) / 45
    return RebarT(b, h, bw, T, dp, ds, de, rec, So, ang)


def REBAREL1(b,  h,  bw,      T,  dp,  ds,      de,  rec,  So,   ang=0):
    ang = -ang * math.atan(1) / 45
    return RebarL1(b, h, bw, T, dp, ds, de, rec, So, ang)


def REBARELE(b,  h,      T,  dp,  ds,      de,  rec,  So,   ang=0):
    ang = -ang * math.atan(1) / 45
    return RebarL(b, h, T, dp, ds, de, rec, So, ang)


def VERTCOMPRES(XY, yx):
    u = RangeToArray(XY)
    return vertCompress(u, yx)


def MOMENTSPUNT(XYA):
    u = RangeToArray(XYA)
    return momentPunt(u)


def MOMENTSPOLY(XY):
    u = RangeToArray(XY)
    return momentPoly(u)


def concrR(b,  h,  ang):
    arry[0][1] = b
    arry[1][2] = h
    arry[0][2] = b
    arry[1][3] = h
    if ang != 0:
        arry = rotatePoint(arry, ang)
    return arry


def concrT2(b,  h,  bw,      T,  b1,  t1,      ang):
    arry[0][1] = b
    arry[1][2] = T
    arry[0][2] = b
    arry[1][3] = T
    arry[0][3] = b1 + bw
    arry[1][4] = h
    arry[0][4] = b1 + bw
    arry[1][5] = h
    arry[0][5] = b1
    arry[1][6] = t1
    arry[0][6] = b1
    arry[1][7] = t1
    if ang != 0:
        arry = rotatePoint(arry, ang)
    return arry


def concrT1(b,  h,  bw,      T,  b1,  ang):
    arry[0][1] = b
    arry[1][2] = T
    arry[0][2] = b
    arry[1][3] = T
    arry[0][3] = b1 + bw
    arry[1][4] = h
    arry[0][4] = b1 + bw
    arry[1][5] = h
    arry[0][5] = b1
    arry[1][6] = T
    arry[0][6] = b1
    arry[1][7] = T
    if ang != 0:
        arry = rotatePoint(arry, ang)
    return arry


def concrT(b,  h,  bw,      T,  ang):
    arry[0][1] = b
    arry[1][2] = T
    arry[0][2] = b
    arry[1][3] = T
    arry[0][3] = (b + bw) / 2
    arry[1][4] = h
    arry[0][4] = (b + bw) / 2
    arry[1][5] = h
    arry[0][5] = (b - bw) / 2
    arry[1][6] = T
    arry[0][6] = (b - bw) / 2
    arry[1][7] = T
    if ang != 0:
        arry = rotatePoint(arry, ang)
    return arry


def concrL1(b,  h,  bw,      T,  ang):
    arry[0][1] = b
    arry[1][2] = T
    arry[0][2] = b
    arry[1][3] = T
    arry[0][3] = bw
    arry[1][4] = h
    arry[0][4] = bw
    arry[1][5] = h
    if ang != 0:
        arry = rotatePoint(arry, ang)
    return arry


def concrL(b,  h,      T,  ang):
    arry[0][1] = b
    arry[1][2] = T
    arry[0][2] = b
    arry[1][3] = T
    arry[0][3] = T
    arry[1][4] = h
    arry[0][4] = T
    arry[1][5] = h
    if ang != 0:
        arry = rotatePoint(arry, ang)
    return arry


def RebarRN(b,  h,  nx,  ny, dp,  ds,  de,  rec,  ang):
    z = Cover(dp, ds, de, rec)
    asp = dp ^ 2 * math.atan(1)
    ass = ds ^ 2 * math.atan(1)
    xyso[0][0] = z[0]
    xyso[1][0] = z[0]
    xyso[2][0] = asp
    xyso[0][1] = b - z[0]
    xyso[1][1] = z[0]
    xyso[2][1] = asp
    xyso[0][2] = z[0]
    xyso[1][2] = h - z[0]
    xyso[2][2] = asp
    xyso[0][3] = b - z[0]
    xyso[1][3] = h - z[0]
    xyso[2][3] = asp
    if nx > 2:
        s1 = (b - 2 * z[2]) / (nx - 1)
        for i in range(1, nx - 2):
            xyso.append([i * s1 + z[2], z[1], ass])
            xyso.append([i * s1 + z[2], h - z[1], ass])
    if ny > 2:
        s2 = (h - 2 * z[2]) / (ny - 1)
        for i in range(1, ny - 2):
            xyso.append([z[1], i * s2 + z[2], ass])
            xyso.append([b - z[1], i * s2 + z[2], ass])
    if ang != 0:
        xyso = rotatePoint(xyso, ang)
    return xyso


def RebarRS(b,  h,      dp,  ds,  de,      rec,  So,  ang):
    z = Cover(dp, ds, de, rec)
    asp = dp ^ 2 * math.atan(1)
    xyso[0][0] = z[0]
    xyso[1][0] = z[0]
    xyso[2][0] = asp
    xyso[0][1] = b - z[0]
    xyso[1][1] = z[0]
    xyso[2][1] = asp
    xyso[0][2] = z[0]
    xyso[1][2] = h - z[0]
    xyso[2][2] = asp
    xyso[0][3] = b - z[0]
    xyso[1][3] = h - z[0]
    xyso[2][3] = asp
    ass = ds ^ 2 * math.atan(1)
    nx = math.ceil((b - 2 * z[2]) / (So + ds) + 1)
    if nx > 2:
        s1 = (b - 2 * z[2]) / (nx - 1)
        for i in range(1, nx - 2):
            xyso.append([i * s1 + z[2], z[1], ass])
            xyso.append([i * s1 + z[2], h - z[1], ass])
    ny = math.ceil((h - 2 * z[2]) / (So + ds) + 1)
    if ny > 2:
        s2 = (h - 2 * z[2]) / (ny - 1)
        for i in range(1, ny - 2):
            xyso.append([z[1], i * s2 + z[2], ass])
            xyso.append([b - z[1], i * s2 + z[2], ass])
    if ang != 0:
        xyso = rotatePoint(xyso, ang)
    return xyso


def RebarT2(b,  h,  bw,      T,  b1,  t1,  dp,      ds,  de,  rec,  So,      ang):
    z = Cover(dp, ds, de, rec)
    asp = dp ^ 2 * math.atan(1)
    ass = ds ^ 2 * math.atan(1)
    xyso[0][0] = z[0]
    xyso[1][0] = z[0]
    xyso[2][0] = asp
    xyso[0][1] = b1 + z[0]
    xyso[1][1] = z[0]
    xyso[2][1] = asp
    xyso[0][2] = b1 + bw - z[0]
    xyso[1][2] = z[0]
    xyso[2][2] = asp
    xyso[0][3] = b - z[0]
    xyso[1][3] = z[0]
    xyso[2][3] = asp
    xyso[0][4] = b - z[0]
    xyso[1][4] = T - z[0]
    xyso[2][4] = asp
    xyso[0][5] = b1 + bw - z[0]
    xyso[1][5] = T - z[0]
    xyso[2][5] = asp
    xyso[0][6] = b1 + bw - z[0]
    xyso[1][6] = h - z[0]
    xyso[2][6] = asp
    xyso[0][7] = b1 + z[0]
    xyso[1][7] = h - z[0]
    xyso[2][7] = asp
    xyso[0][8] = b1 + z[0]
    xyso[1][8] = t1 - z[0]
    xyso[2][8] = asp
    xyso[0][9] = z[0]
    xyso[1][9] = t1 - z[0]
    xyso[2][9] = asp
    w = []
    r = []
    q = []
    pos = []
    w[0] = b1 + z[1] - z[2]
    r[0] = z[2]
    q[0] = z[1]
    pos[0] = True
    w[1] = bw - 2 * z[2]
    r[1] = b1 + z[2]
    q[1] = z[1]
    pos[1] = True
    w[2] = b - b1 - bw + z[1] - z[2]
    r[2] = b1 + bw - z[1]
    q[2] = z[1]
    pos[2] = True
    w[3] = T - 2 * z[2]
    r[3] = z[2]
    q[3] = b - z[1]
    pos[3] = False
    w[4] = b - b1 - bw + z[1] - z[2]
    r[4] = b1 + bw - z[1]
    q[4] = T - z[1]
    pos[4] = True
    w[5] = h - T + z[1] - z[2]
    r[5] = T - z[1]
    q[5] = b1 + bw - z[1]
    pos[5] = False
    w[6] = bw - 2 * z[2]
    r[6] = b1 + z[2]
    q[6] = h - z[1]
    pos[6] = True
    w[7] = h - t1 + z[1] - z[2]
    r[7] = t1 - z[1]
    q[7] = b1 + z[1]
    pos[7] = False
    w[8] = b1 + z[1] - z[2]
    r[8] = z[2]
    q[8] = t1 - z[1]
    pos[8] = True
    w[9] = t1 - 2 * z[2]
    r[9] = z[2]
    q[9] = z[1]
    pos[9] = False
    n = []
    s = []
    for i in range(9):
        n[i] = math.ceil(w[i] / (So + ds) - 1)
        s[i] = w[i] / (n[i] + 1)
    for j in range(9):
        for i in range(1, n[j]):
            if pos(j):
                xyso.append([i * s(j) + r(j), q(j), ass])
            else:
                xyso.append([q(j), i * s(j) + r(j), ass])
    if ang != 0:
        xyso = rotatePoint(xyso, ang)
    return xyso


def RebarT1(b,  h,  bw,      T,  b1,  dp,      ds,  de,  rec,  So,      ang):
    z = Cover(dp, ds, de, rec)
    asp = dp ^ 2 * math.atan[1]
    ass = ds ^ 2 * math.atan[1]
    xyso[0][0] = z[0]
    xyso[1][0] = z[0]
    xyso[2][0] = asp
    xyso[0][1] = b1 + z[0]
    xyso[1][1] = z[0]
    xyso[2][1] = asp
    xyso[0][2] = b1 + bw - z[0]
    xyso[1][2] = z[0]
    xyso[2][2] = asp
    xyso[0][3] = b - z[0]
    xyso[1][3] = z[0]
    xyso[2][3] = asp
    xyso[0][4] = b - z[0]
    xyso[1][4] = T - z[0]
    xyso[2][4] = asp
    xyso[0][5] = b1 + bw - z[0]
    xyso[1][5] = T - z[0]
    xyso[2][5] = asp
    xyso[0][6] = b1 + bw - z[0]
    xyso[1][6] = h - z[0]
    xyso[2][6] = asp
    xyso[0][7] = b1 + z[0]
    xyso[1][7] = h - z[0]
    xyso[2][7] = asp
    xyso[0][8] = b1 + z[0]
    xyso[1][8] = T - z[0]
    xyso[2][8] = asp
    xyso[0][9] = z[0]
    xyso[1][9] = T - z[0]
    xyso[2][9] = asp
    w = []
    r = []
    q = []
    pos = []
    w[0] = b1 + z[1] - z[2]
    r[0] = z[2]
    q[0] = z[1]
    pos[0] = True
    w[1] = bw - 2 * z[2]
    r[1] = b1 + z[2]
    q[1] = z[1]
    pos[1] = True
    w[2] = b - b1 - bw + z[1] - z[2]
    r[2] = b1 + bw - z[1]
    q[2] = z[1]
    pos[2] = True
    w[3] = T - 2 * z[2]
    r[3] = z[2]
    q[3] = b - z[1]
    pos[3] = False
    w[4] = b - b1 - bw + z[1] - z[2]
    r[4] = b1 + bw - z[1]
    q[4] = T - z[1]
    pos[4] = True
    w[5] = h - T + z[1] - z[2]
    r[5] = T - z[1]
    q[5] = b1 + bw - z[1]
    pos[5] = False
    w[6] = bw - 2 * z[2]
    r[6] = b1 + z[2]
    q[6] = h - z[1]
    pos[6] = True
    w[7] = h - T + z[1] - z[2]
    r[7] = T - z[1]
    q[7] = b1 + z[1]
    pos[7] = False
    w[8] = b1 + z[1] - z[2]
    r[8] = z[2]
    q[8] = T - z[1]
    pos[8] = True
    w[9] = T - 2 * z[2]
    r[9] = z[2]
    q[9] = z[1]
    pos[9] = False
    n = []
    s = []
    for i in range[9]:
        n[i] = math.ceil(w[i] / (So + ds) - 1)
        s[i] = w[i] / (n[i] + 1)
    for j in range[9]:
        for i in range(1, n[j]):
            if pos(j):
                xyso.append([i * s(j) + r(j), q(j), ass])
            else:
                xyso.append([q(j), i * s(j) + r(j), ass])
    if ang != 0:
        xyso = rotatePoint(xyso, ang)
    return xyso


def RebarT(b,  h,  bw,      T,  dp,  ds,  de,      rec,  So,  ang):
    z = Cover(dp, ds, de, rec)
    asp = dp ^ 2 * math.atan[1]
    ass = ds ^ 2 * math.atan[1]
    xyso[0][0] = z[0]
    xyso[1][0] = z[0]
    xyso[2][0] = asp
    xyso[0][1] = (b - bw) / 2 + z[0]
    xyso[1][1] = z[0]
    xyso[2][1] = asp
    xyso[0][2] = (b + bw) / 2 - z[0]
    xyso[1][2] = z[0]
    xyso[2][2] = asp
    xyso[0][3] = b - z[0]
    xyso[1][3] = z[0]
    xyso[2][3] = asp
    xyso[0][4] = b - z[0]
    xyso[1][4] = T - z[0]
    xyso[2][4] = asp
    xyso[0][5] = (b + bw) / 2 - z[0]
    xyso[1][5] = T - z[0]
    xyso[2][5] = asp
    xyso[0][6] = (b + bw) / 2 - z[0]
    xyso[1][6] = h - z[0]
    xyso[2][6] = asp
    xyso[0][7] = (b - bw) / 2 + z[0]
    xyso[1][7] = h - z[0]
    xyso[2][7] = asp
    xyso[0][8] = (b - bw) / 2 + z[0]
    xyso[1][8] = T - z[0]
    xyso[2][8] = asp
    xyso[0][9] = z[0]
    xyso[1][9] = T - z[0]
    xyso[2][9] = asp
    w = []
    r = []
    q = []
    pos = []
    w[0] = (b - bw) / 2 + z[1] - z[2]
    r[0] = z[2]
    q[0] = z[1]
    pos[0] = True
    w[1] = bw - 2 * z[2]
    r[1] = (b - bw) / 2 + z[2]
    q[1] = z[1]
    pos[1] = True
    w[2] = (b - bw) / 2 + z[1] - z[2]
    r[2] = (b + bw) / 2 - z[1]
    q[2] = z[1]
    pos[2] = True
    w[3] = T - 2 * z[2]
    r[3] = z[2]
    q[3] = b - z[1]
    pos[3] = False
    w[4] = (b - bw) / 2 + z[1] - z[2]
    r[4] = (b + bw) / 2 - z[1]
    q[4] = T - z[1]
    pos[4] = True
    w[5] = h - T + z[1] - z[2]
    r[5] = T - z[1]
    q[5] = (b + bw) / 2 - z[1]
    pos[5] = False
    w[6] = bw - 2 * z[2]
    r[6] = (b - bw) / 2 + z[2]
    q[6] = h - z[1]
    pos[6] = True
    w[7] = h - T + z[1] - z[2]
    r[7] = T - z[1]
    q[7] = (b - bw) / 2 + z[1]
    pos[7] = False
    w[8] = (b - bw) / 2 + z[1] - z[2]
    r[8] = z[2]
    q[8] = T - z[1]
    pos[8] = True
    w[9] = T - 2 * z[2]
    r[9] = z[2]
    q[9] = z[1]
    pos[9] = False
    n = []
    s = []
    for i in range[9]:
        n[i] = math.ceil(w[i] / (So + ds) - 1)
        s[i] = w[i] / (n[i] + 1)
    for j in range[9]:
        for i in range(1, n[j]):
            if pos(j):
                xyso.append([i * s(j) + r(j), q(j), ass])
            else:
                xyso.append([q(j), i * s(j) + r(j), ass])
    if ang != 0:
        xyso = rotatePoint(xyso, ang)
    return xyso


def RebarL1(b,  h,  bw,      T,  dp,      ds,  de,  rec,  So,      ang):
    z = Cover(dp, ds, de, rec)
    asp = dp ^ 2 * math.math.atan(1)
    ass = ds ^ 2 * math.math.atan(1)
    xyso[0][0] = z[0]
    xyso[1][0] = z[0]
    xyso[2][0] = asp
    xyso[0][1] = bw - z[0]
    xyso[1][1] = z[0]
    xyso[2][1] = asp
    xyso[0][2] = b - z[0]
    xyso[1][2] = z[0]
    xyso[2][2] = asp
    xyso[0][3] = b - z[0]
    xyso[1][3] = T - z[0]
    xyso[2][3] = asp
    xyso[0][4] = bw - z[0]
    xyso[1][4] = T - z[0]
    xyso[2][4] = asp
    xyso[0][5] = bw - z[0]
    xyso[1][5] = h - z[0]
    xyso[2][5] = asp
    xyso[0][6] = z[0]
    xyso[1][6] = h - z[0]
    xyso[2][6] = asp
    xyso[0][7] = z[0]
    xyso[1][7] = T - z[0]
    xyso[2][7] = asp
    w = []
    r = []
    q = []
    pos = []
    w[0] = bw - 2 * z[2]
    r[0] = z[2]
    q[0] = z[1]
    pos[0] = True
    w[1] = b - bw + z[1] - z[2]
    r[1] = bw - z[1]
    q[1] = z[1]
    pos[1] = True
    w[2] = T - 2 * z[2]
    r[2] = z[2]
    q[2] = b - z[1]
    pos[2] = False
    w[3] = b - bw + z[1] - z[2]
    r[3] = bw - z[1]
    q[3] = T - z[1]
    pos[3] = True
    w[4] = h - T + z[1] - z[2]
    r[4] = T - z[1]
    q[4] = bw - z[1]
    pos[4] = False
    w[5] = bw - 2 * z[2]
    r[5] = z[2]
    q[5] = h - z[1]
    pos[5] = True
    w[6] = h - T + z[1] - z[2]
    r[6] = T - z[1]
    q[6] = z[1]
    pos[6] = False
    w[7] = T - 2 * z[2]
    r[7] = z[2]
    q[7] = z[1]
    pos[7] = False
    n = []
    s = []
    for i in range(7):
        n[i] = math.ceil(w[i] / (So + ds) - 1)
        s[i] = w[i] / (n[i] + 1)
    for j in range(7):
        for i in range(1, n[j]):
            if pos(j):
                xyso.append([i * s(j) + r(j), q(j), ass])
            else:
                xyso.append([q(j), i * s(j) + r(j), ass])
    if ang != 0:
        xyso = rotatePoint(xyso, ang)
    return xyso


def RebarL(b,  h,      T,  dp,      ds,  de,  rec,  So,      ang):
    z = Cover(dp, ds, de, rec)
    asp = dp ^ 2 * math.atan[1]
    ass = ds ^ 2 * math.atan[1]
    xyso = []
    xyso[0][0] = z[0]
    xyso[1][0] = z[0]
    xyso[2][0] = asp
    xyso[0][1] = T - z[0]
    xyso[1][1] = z[0]
    xyso[2][1] = asp
    xyso[0][2] = b - z[0]
    xyso[1][2] = z[0]
    xyso[2][2] = asp
    xyso[0][3] = b - z[0]
    xyso[1][3] = T - z[0]
    xyso[2][3] = asp
    xyso[0][4] = T - z[0]
    xyso[1][4] = T - z[0]
    xyso[2][4] = asp
    xyso[0][5] = T - z[0]
    xyso[1][5] = h - z[0]
    xyso[2][5] = asp
    xyso[0][6] = z[0]
    xyso[1][6] = h - z[0]
    xyso[2][6] = asp
    xyso[0][7] = z[0]
    xyso[1][7] = T - z[0]
    xyso[2][7] = asp
    w = []
    r = []
    q = []
    pos = []
    w[0] = T - 2 * z[2]
    r[0] = z[2]
    q[0] = z[1]
    pos[0] = True
    w[1] = b - T + z[1] - z[2]
    r[1] = T - z[1]
    q[1] = z[1]
    pos[1] = True
    w[2] = T - 2 * z[2]
    r[2] = z[2]
    q[2] = b - z[1]
    pos[2] = False
    w[3] = b - T + z[1] - z[2]
    r[3] = T - z[1]
    q[3] = T - z[1]
    pos[3] = True
    w[4] = h - T + z[1] - z[2]
    r[4] = T - z[1]
    q[4] = T - z[1]
    pos[4] = False
    w[5] = T - 2 * z[2]
    r[5] = z[2]
    q[5] = h - z[1]
    pos[5] = True
    w[6] = h - T + z[1] - z[2]
    r[6] = T - z[1]
    q[6] = z[1]
    pos[6] = False
    w[7] = T - 2 * z[2]
    r[7] = z[2]
    q[7] = z[1]
    pos[7] = False
    n = []
    s = []
    for i in range(7):
        n[i] = math.ceil(w[i] / (So + ds) - 1)
        s[i] = w[i] / (n[i] + 1)
    for j in range(7):
        for i in range(1, n[j]):
            if pos(j):
                xyso.append([i * s(j) + r(j), q(j), ass])
            else:
                xyso.append([q(j), i * s(j) + r(j), ass])
    if ang != 0:
        xyso = rotatePoint(xyso, ang)
    return xyso


def diagraFlexCompress(xygo, xyso, alfa, fc, fy,  Es,  fcp):
    sc = -fatiga * fc
    ey = fy / Es * 1000
    fyc = fy / sc
    cplast = centroidPlast(xygo, xyso, fyc, fy, sc)
    XYG = trasltPlastg(xygo, cplast)
    XYS = trasltPlasts(xyso, cplast)
    xyga = rotatePoint(XYG, alfa)
    xysa = rotatePoint(XYS, alfa)
    ycomp = max(xyga[2])
    cmax = ycomp - min(xyga[2])
    ks = len(xysa[2])
    d = []
    for i in range(ks):
        d[i] = ycomp - xysa[1][i]
    dtens = max(d)
    dcomp = min(d)
    xt = dtens / (1 - ety / ecu)
    etc = ety / ecu
    eyc = ey / ecu
    cd = []
    cd[0] = (ey * dtens - ety * dcomp) / (ey - ety)
    cd[1] = 0
    cd[2] = 0.5 * xt
    cd[3] = 0.75 * xt
    cd[4] = xt
    cd[5] = dtens / (1 - 0.75 * etc)
    cd[6] = dtens / (1 - 0.5 * etc)
    cd[7] = dtens / (1 - eyc)
    cd[8] = dtens / (1 - 0.75 * eyc)
    cd[9] = dtens / (1 - 0.5 * eyc)
    cd[10] = dtens / (1 - 0.25 * eyc)
    cd[11] = dtens
    cd[12] = dtens / 0.9
    cd[13] = dtens / (1 + eyc)
    arry = []
    for i in range(13):
        arry[i] = flector(xyga, xysa, d, dtens, cmax, ycomp,
                          alfa, cd[i], fcp, cplast[2], sc, ey, Es, fy, fc)
    return arry


def Cover(dp,  ds,  de,      rec):
    arry = []
    arry[0] = rec + de + dp / 2
    arry[1] = rec + de + ds / 2
    arry[2] = arry[0] + (dp - ds) / 2
    return arry


def xCompressGral(xygo,  xyso,  alfa, fc,  fy,  Es,  P, fcp):
    sc = -fatiga * fc
    ey = fy / Es * 1000
    fyc = fy / sc
    cplast = centroidPlast(xygo, xyso, fyc, fy, sc)
    XYG = trasltPlastg(xygo, cplast)
    XYS = trasltPlasts(xyso, cplast)
    xyga = rotatePoint(XYG, alfa)
    xysa = rotatePoint(XYS, alfa)
    ycomp = max(xyga[2])
    cmax = ycomp - min(xyga[2])
    u1 = -cmax
    u2 = 2 * cmax
    for i in range(15):
        if abs(u1 - u2) < 0.0001:
            break
        c = (u1 + u2) / 2
        f1 = P - Normal(xyga, xysa, cmax, ycomp, c, fcp, sc, ey, Es, fy, fc)
        if f1 < 0:
            u1 = c
        else:
            u2 = c
    return c


def xCompressGralsu(xygo,  xyso,  alfa, fc,  fy,  Es,  P, fcp):
    sc = -fatiga * fc
    ey = fy / Es * 1000
    fyc = fy / sc
    cplast = centroidPlast(xygo, xyso, fyc, fy, sc)
    XYG = trasltPlastg(xygo, cplast)
    XYS = trasltPlasts(xyso, cplast)
    xyga = rotatePoint(XYG, alfa)
    xysa = rotatePoint(XYS, alfa)
    ycomp = max(xyga[2])
    cmax = ycomp - min(xyga[2])
    u1 = -cmax
    u2 = 2 * cmax
    P = max(P, cplast[2])
    for i in range(15):
        if abs(u1 - u2) < 0.0001:
            break
        c = (u1 + u2) / 2
        f1 = P - Normalsu(xyga, xysa, cmax, ycomp, c, fcp,
                          cplast[2], sc, ey, Es, fy, fc)
        if f1 < 0:
            u1 = c
        else:
            u2 = c
    return c


def FlexCompressGral(xygo,  xyso,  alfa,      c,  fc,  fy,  Es,      fcp):
    sc = -fatiga * fc
    ey = fy / Es * 1000
    fyc = fy / sc
    cplast = centroidPlast(xygo, xyso, fyc, fy, sc)
    XYG = trasltPlastg(xygo, cplast)
    XYS = trasltPlasts(xyso, cplast)
    xyga = rotatePoint(XYG, alfa)
    xysa = rotatePoint(XYS, alfa)
    ycomp = max(xyga[2])
    cmax = ycomp - min(xyga[2])
    ks = len(xysa[2])
    d = []
    for i in range(ks):
        d[i] = ycomp - xysa[1][i]
    dtens = max(d)
    return flector(xyga, xysa, d, dtens, cmax, ycomp, alfa, c, fcp, cplast[2], sc, ey, Es, fy, fc)


def FlexCompressGralsu(xygo,  xyso,  alfa,      c,  fc,  fy,  Es,      fcp):
    sc = -fatiga * fc
    ey = fy / Es * 1000
    fyc = fy / sc
    cplast = centroidPlast(xygo, xyso, fyc, fy, sc)
    XYG = trasltPlastg(xygo, cplast)
    XYS = trasltPlasts(xyso, cplast)
    xyga = rotatePoint(XYG, alfa)
    xysa = rotatePoint(XYS, alfa)
    ycomp = max(xyga[2])
    cmax = ycomp - min(xyga[2])
    ks = len(xysa[2])
    d = []
    for i in range(ks):
        d[i] = ycomp - xysa[1][i]
    dtens = max(d)
    return flectorsu(xyga, xysa, d, dtens, cmax, ycomp, alfa, c, fcp, sc, ey, Es, fy, fc)


def Normal(xyga, xysa, cmax,  ycomp,  c,  fcp,  sc,  ey,  Es,  fy,  fc, ks):
    d = []
    for i in range(ks):
        d.append(ycomp - xysa[1][i])
    dtens = max(d)
    ess = []
    if c < dtens / (1 - ety / ecu):
        for i in range(ks):
            ess.append((c - d[i]) / (c - dtens) * ety)
    else:
        for i in range(ks):
            ess.append((1 - d[i] / c) * ecu)
    esmax = max(ess)
    if esmax <= ey:
        phiu = phic
    else:
        if esmax >= ety:
            phiu = phit
        else:
            phiu = phic + (phit - phic) * (esmax - ey) / (ety - ey)
    if c > 0:
        if fc < fcp:
            a = 0.85 * c
        else:
            if fc > 2 * fcp:
                a = 0.65 * c
            else:
                a = (0.85 - (0.85 - 0.65) * (fc / fcp - 1)) * c
        a = min(a, cmax)
    else:
        a = 0
    ss = []
    nso = []
    for i in range(ks):
        if abs(ess[i]) < ey:
            ss.append(ess[i] * Es / 1000)
        else:
            ss.append(math.copysign(fy, ess[i]))
        if d[i] < a:
            ss[i] = ss[i] - sc
        nso.append(ss[i] * xysa[2][i] * phiu)
    ns = sum(nso)
    if c > 0:
        xyc = xyga
        if a < cmax:
            xyc = vertCompress(xyc, ycomp - a)
        ng = areaPoly(xyc) * sc * phiu
    return ns + ng


def flector(xyga,  xysa,  d,  dtens,  cmax,      ycomp,  alfa,  c,  fcp,  nmax,      sc,  ey,  Es,  fy,  fc):
    ks = len(xysa[2])
    ess = []
    if c < dtens / (1 - ety / ecu):
        for i in range(ks):
            ess.append((c - d[i]) / (c - dtens) * ety)
    else:
        for i in range(ks):
            ess.append((1 - d[i] / c) * ecu)
    esmax = max(ess)
    if esmax <= ey:
        phiu = phic
    else:
        if esmax >= ety:
            phiu = phit
        else:
            phiu = phic + (phit - phic) * (esmax - ey) / (ety - ey)
    if c > 0:
        if fc < fcp:
            a = 0.85 * c
        else:
            if fc > 2 * fcp:
                a = 0.65 * c
            else:
                a = (0.85 - (0.85 - 0.65) * (fc / fcp - 1)) * c
        a = min(a, cmax)
    else:
        a = 0
    ss = []
    nso = []
    msxya = []
    nmsxy = []
    for i in range(ks):
        if abs(ess[i]) < ey:
            ss[i] = ess[i] * Es / 1000
        else:
            ss[i] = math.copysign(ess[i], fy)
        if d[i] < a:
            ss[i] = ss[i] - sc
        nso[i] = ss[i] * xysa[2][i] * phiu
        msxya[0][i] = nso[i] * xysa[1][i]
        msxya[1][i] = nso[i] * xysa[0][i]
        nmsxy[0] = sum(nso)
        nmsxy[1] = sum(msxya[1])
        nmsxy[2] = sum(msxya[2])
    if c > 0:
        xyc = xyga
        if a < cmax:
            xyc = vertCompress(xyc, ycomp - a)
        nmgxy = momentPoly(xyc)
        nmgxy[0] = nmgxy[0] * sc * phiu
        nmgxy[1] = nmgxy[1] * sc * phiu
        nmgxy[2] = nmgxy[2] * sc * phiu
    Mnxa = nmgxy[1] + nmsxy[1]
    Mnya = nmgxy[2] + nmsxy[2]
    arry = []
    arry[0] = max(nmgxy[0] + nmsxy[0], nmax)
    arry[1] = -Mnxa * math.sin(alfa) + Mnya * math.cos(alfa)
    arry[2] = -Mnxa * math.cos(alfa) - Mnya * math.sin(alfa)
    return arry


def flectorsu(xyga,  xysa,  d,  dtens,  cmax,      ycomp,  alfa,  c,  fcp,      sc,  ey,  Es,  fy,  fc):
    ks = len(xysa[2])
    ess = []
    if c < dtens / (1 - ety / ecu):
        for i in range(ks):
            ess[i] = (c - d[i]) / (c - dtens) * ety
    else:
        for i in range(ks):
            ess[i] = (1 - d[i] / c) * ecu
    esmax = max(ess)
    if esmax <= ey:
        phiu = phic
    else:
        if esmax >= ety:
            phiu = phit
        else:
            phiu = phic + (phit - phic) * (esmax - ey) / (ety - ey)
    if c > 0:
        if fc < fcp:
            a = 0.85 * c
        else:
            if fc > 2 * fcp:
                a = 0.65 * c
            else:
                a = (0.85 - (0.85 - 0.65) * (fc / fcp - 1)) * c
        a = min(a, cmax)
    else:
        a = 0
    ss = []
    nso = []
    msxya = []
    nmsxy = []
    for i in range(ks):
        if abs(ess[i]) < ey:
            ss[i] = ess[i] * Es / 1000
        else:
            ss[i] = math.copysign(ess[i], fy)
        if d[i] < a:
            ss[i] = ss[i] - sc
        nso[i] = ss[i] * xysa[2][i] * phiu
        msxya[0][i] = nso[i] * xysa[1][i]
        msxya[1][i] = nso[i] * xysa[0][i]
        nmsxy[0] = sum(nso)
        nmsxy[1] = sum(msxya[1])
        nmsxy[2] = sum(msxya[2])
    if c > 0:
        xyc = xyga
        if a < cmax:
            xyc = vertCompress(xyc, ycomp - a)
        nmgxy = momentPoly(xyc)
        nmgxy[0] = nmgxy[0] * sc * phiu
        nmgxy[1] = nmgxy[1] * sc * phiu
        nmgxy[2] = nmgxy[2] * sc * phiu
    Mnxa = nmgxy[1] + nmsxy[1]
    Mnya = nmgxy[2] + nmsxy[2]
    arry = []
    arry[0] = nmgxy[0] + nmsxy[0]
    arry[1] = -Mnxa * math.sin(alfa) + Mnya * math.cos(alfa)
    arry[2] = -Mnxa * math.cos(alfa) - Mnya * math.sin(alfa)
    return arry


def Normalsu(xyga,  xysa,  cmax, ycomp,  c,  fcp,  nmax,      sc,  ey,  Es,  fy,  fc):
    ks = len(xysa[2])
    d = []
    for i in range(ks):
        d.append(ycomp - xysa[1][i])
    dtens = max(d)
    ess = []
    if c < dtens / (1 - ety / ecu):
        for i in range(ks):
            ess[i] = (c - d[i]) / (c - dtens) * ety
    else:
        for i in range(ks):
            ess[i] = (1 - d[i] / c) * ecu
    esmax = max(ess)
    if esmax <= ey:
        phiu = phic
    else:
        if esmax >= ety:
            phiu = phit
        else:
            phiu = phic + (phit - phic) * (esmax - ey) / (ety - ey)
    if c > 0:
        if fc < fcp:
            a = 0.85 * c
        else:
            if fc > 2 * fcp:
                a = 0.65 * c
            else:
                a = (0.85 - (0.85 - 0.65) * (fc / fcp - 1)) * c
        a = min(a, cmax)
    else:
        a = 0
    ss = []
    nso = []
    for i in range(ks):
        if abs(ess[i]) < ey:
            ss[i] = ess[i] * Es / 1000
        else:
            ss[i] = math.copysign(ess[i], fy)
        if d[i] < a:
            ss[i] = ss[i] - sc
        nso[i] = ss[i] * xysa[2][i] * phiu
    ns = sum(nso)
    if c > 0:
        xyc = xyga
        if a < cmax:
            xyc = vertCompress(xyc, ycomp - a)
        ng = areaPoly(xyc) * sc * phiu
    return max(ns + ng, nmax)


def nvig(sh, a,  k):
    # sh espaciamiento entre barras
    # a=sqrt(wt*fy/fc*b*h/pi)
    # k=b+sh-2*rec-2*de
    u1 = 2
    for iloop in range[9]:
        u2 = u1 - (a * u1 + sh * u1 ^ 2 - k) / (a + 2 * sh * u1)
        if u1 == u2:
            break
        u1 = u2
    nvig = math.ceil(u2 ^ 2)


def centroidPlast(xygo, xyso, fyc, fy, sc):
    arry = [0, 0, 0]
    mg = momentPoly(xygo)
    ms = momentPunt(xyso)
    arry[0] = (mg[2] - (fyc + 1) * ms[2]) / (mg[0] - (fyc + 1) * ms[0])
    arry[1] = (mg[1] - (fyc + 1) * ms[1]) / (mg[0] - (fyc + 1) * ms[0])
    arry[2] = 0.8 * phic * (sc * (mg[0] - ms[0]) - fy * ms[0])
    return arry


def trasltPlastg(xygo,  cplas, k):
    arry = [[], []]
    for i in range(k):
        arry[0].append(xygo[0][i] - cplas[0])
        arry[1].append(xygo[1][i] - cplas[1])
    return arry


def trasltPlasts(xyso,  cplas, k):
    arry = xyso
    for i in range(k):
        arry[0][i] = xyso[0][i] - cplas[0]
        arry[1][i] = xyso[1][i] - cplas[1]
    return arry


def vertCompress(xyc,  yx):
    kg = len(xyc[2])
    for i in range(kg - 1, 0, -1):
        if ((yx > xyc[1][i]) and (yx < xyc[1][i + 1])) or ((yx < xyc[1][i]) and (yx > xyc[1][i + 1])):
            xyc.insert[i + 1][xyc[0][i + 1] + (xyc[0][i] - xyc[0][i + 1]) / (
                xyc[1][i] - xyc[1][i + 1]) * (yx - xyc[1][i + 1]), yx]
    for i in range(len(xyc[2]), 0, -1):
        if xyc[1][i] < yx:
            xyc.pop(i)
    xyc[0].append(xyc[0][0])
    xyc[1].append(xyc[1][0])
    return xyc


def momentPoly(xygo, kg):
    arry = [0, 0, 0]
    asi = []
    for i in range(kg):
        asi.append(xygo[0][i] * xygo[1][i + 1] - xygo[0][i + 1] * xygo[1][i])
        arry[0] += asi[i]
        arry[1] += (xygo[1][i] + xygo[1][i + 1]) * asi[i]
        arry[2] += (xygo[0][i] + xygo[0][i + 1]) * asi[i]
    return [arry[0] / 2, arry[1] / 6, arry[2] / 6]


def areaPoly(xygo, k):
    areaPoly = 0
    for i in range(k):
        areaPoly += xygo[0][i] * xygo[1][i + 1] - xygo[0][i + 1] * xygo[1][i]
    return areaPoly / 2


def momentPunt(xyso, k):
    arry = [0, 0, 0]
    for i in range(k):
        arry[0] += xyso[2][i]
        arry[1] += xyso[1][i] * xyso[2][i]
        arry[2] += xyso[0][i] * xyso[2][i]
    return arry


def rotatePoint(xy, z):
    xya = xy
    cosz = math.cos(z)
    sinz = math.sin(z)
    for i in range(len(xya[0])):
        xya[0][i] = xy[0][i] * cosz + xy[1][i] * sinz
        xya[1][i] = -xy[0][i] * sinz + xy[1][i] * cosz
    return xya


def RangeToArray(XY):
    n = XY.Columns.Count - 1
    m = XY.Rows.Count - 1
    arry = []
    for i in range(n):
        for j in range(m):
            arry[i][j] = XY(j + 1, i + 1)
    return arry


def main():
    kg = 4
    xygo = [[0, 9, 26, 7, 0], [0, -3, 19, 24, 0]]
    ks = 3
    xyso = [[7, 18, 10], [3, 16, 18], [3.142, 2.011, 2.011]]
    print(momentPunt(xyso, ks))
    alfa = 3
    print(rotatePoint(xygo, alfa))
    print(areaPoly(xygo, kg))
    print(momentPoly(xygo, kg))


if __name__ == "__main__":
    main()
