function Type=ObjectDefine_Passive(ObjName);
%This is a function to define the object type according to the object name:
%Last digit: <4, Good(Type1); <8, Bad(Type2), <16;Novel(Type3)
%by Xuefei Yu 08102020
Type=NaN*ones(length(ObjName),1);
Type(ObjName>9999)=3;
type_tmp=Type(isnan(Type));
ObjShort=ObjName(ObjName<=9999);
LastDigit=mod(ObjShort,floor(ObjShort./(ObjShort/10)));

type_tmp(LastDigit>=8)=3;

type_tmp(LastDigit<4)=1;
type_tmp((LastDigit>=4)&(LastDigit<8))=2;

Type(isnan(Type))=type_tmp;

end