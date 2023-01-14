
function m_s_f = func_vd_msf (y)

clear m_s_f;

[B,A] = butter(9,.33,'low');  %.5 or .33?
y1 = filter(B,A,y);

m_s_f=sum(abs(y1));


