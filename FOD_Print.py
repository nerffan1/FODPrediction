from FOD import *
from FFOD import *
from BFOD import *

def PrintSidebySide(a, b):
    print(50*'-')
    print(f"{a.mBold.mName}-{a.mMeek.mName} {type(a)}")
    print("{:<10} {:<13} {:<5} {:<13}".format('', 'Original FOD', '', 'Associated FOD'))
    print(f"{'BoldR:':<10} {a.mBoldR:<13.6f} {'':<5} {b.mBoldR:<13.6f}")
    print(f"{'MeekR:':<10} {a.mMeekR:<13.6f} {'':<5} {b.mMeekR:<13.6f}")
    print(f"{'BoldA:':<10} {np.rad2deg(a.mBoldAngle):<13.6f} {'':<5} {np.rad2deg(b.mBoldAngle):<13.6f}")
    print(f"{'MeekA:':<10} {np.rad2deg(a.mMeekAngle):<13.6f} {'':<5} {np.rad2deg(b.mMeekAngle):<13.6f}")
    print(f"{'Port.:':<10} {a.mBoldPortion:<13.6f} {'':<5} {b.mBoldPortion:<13.6f}")
    print(f"Dist: {dist(a.mPos,b.mPos)}")