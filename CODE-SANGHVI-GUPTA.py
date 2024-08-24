import numpy as np
import matplotlib.pyplot as plt
import math

# Submitted by:
# SANGHVI Mahek Atul & GUPTA Vivek Anantkumar

plt.close('all')

#AIRCRAFT DATA
#The sircraft considered for this project is AIRBUS A320
#all units are in SI units

b = 35.8
rho = 0.310828
u_inf = 248
n = 51
a0 = 6.28
s = 122.6
AR = pow(b,2)/s
weight = 765180

# chord distribution
wingtaper = 0.266
meanchord = s/b
lambdaa = (1-wingtaper)/(1+wingtaper)
chord = np.zeros(n)
y = np.linspace(-b/2,b/2,n)
x = np.linspace(-b/2,b/2,n)
theta = np.linspace(0.000000001,179.9999999999,n)

eta = np.cos(np.radians(theta))
CT = meanchord * ((1+lambdaa)-(2*lambdaa*(np.abs(np.cos(np.radians(theta))))))

N = 2*np.linspace(1,n,n)-1
MATA = np.zeros((n, n))
MATB = np.zeros((n, 1))
MATBncoef = np.zeros((n, 1))

Z = (a0/(4*b))

for i in range(n):
    for j in range(n):
        MATA[i, j] = np.sin((N[j]*np.radians(theta[i]))) + (Z * N[j] * CT[i] * (np.sin(N[j]*np.radians(theta[i])))/(np.sin(np.radians(theta[i]))))

for i in range(n):
    MATB[i] = 2*Z*CT[i]

MATBncoef = np.dot(np.linalg.inv(MATA), MATB)

L = weight # for cruise condtion lift = Weight.
Cl = L / (0.5*rho*(u_inf**2)*s)
Cl = round(Cl,3)

print("Lift coef: ", Cl)


#lift Slope
lift_slope = a0/(1+(a0/(np.pi*AR)))
print("lift_slope: ",lift_slope)
alpha_wing = Cl/lift_slope
print("alpha_wing: ",np.rad2deg(alpha_wing))

#circulationulation
circulation = np.zeros(n)
SUM = np.zeros(n)

for i in range(n):
    for j in range(n):
       SUM[i] = SUM[i] + MATBncoef[j]*np.sin(N[j]*np.radians(theta[i]))
    # print(np.sin(N[j]*np.radians(theta[i])))
    circulation[i] = u_inf * b * (alpha_wing) * SUM[i]


#span loading
SPload = np.zeros(n)
SUM1 = np.zeros(n)
for i in range(n):
    for j in range(n):
       SUM1[i] = SUM1[i] + MATBncoef[j]*np.sin((N[j]*np.radians(theta[i])))
    SPload[i] = 2 * AR * (alpha_wing) * SUM1[i]

#calculation of Downwash_angle angle of attack
Downwash_angle = np.zeros(n)
SUM2 = np.zeros(n)
for i in range(n):
    for j in range(n):
            SUM2[i] = SUM2[i] + 0.5 * (N[j] * MATBncoef[j] * np.sin(N[j] * (np.radians(theta[i]))) / (np.sin((np.radians(theta[i])))))
    Downwash_angle[i] =  (alpha_wing)* SUM2[i]

# calculation of alpha effective
alpha_effe = np.zeros(n)

SUM3 = np.zeros(n)

for i in range(n):

    alpha_effe[i] = ((alpha_wing)) - Downwash_angle[i]
################# 5th question
CDinduced = 0
SUM5 = 0
for j in range(n):
    SUM5 = SUM5 + (N[j] * ((MATBncoef[j])/MATBncoef[0])**2)   
CDinduced = (((np.pi/4)*AR) * (Cl/((np.pi/2)*AR))**2) * SUM5
print("ID is: " , CDinduced)

#####Oswald efficiency
oswald = (Cl**2)/(np.pi*AR* CDinduced)
print("Efficiency oswald: ", oswald)

#########-------> Question 6
circulation_zero = circulation[25]
addition = 0
d_theta = (np.radians(theta[2])-(np.radians(theta[1])))
for i in range(n):
    addition = addition + circulation[i] * b/2 * np.sin(np.radians(theta[i])) * d_theta

I = addition/2
y0 = I/(circulation_zero)
b0 = 2*y0
S = b0/b
print("s:", S)

E0 = (S**2/oswald)*(2/np.pi)*(circulation_zero**2)
power = (((2*np.pi*E0)/(circulation_zero**2)) + 0.5)
e_to_the_x = math.exp(power)
rc = b0/ e_to_the_x
print("rc/b: ",rc/b)
print("circulation zero/U*b", circulation_zero/(u_inf*b))

############-------> Question 7
m = 51 # discretization points

r_right = np.linspace(((-3*b0)/2),(b0/2),m)
r_left = np.linspace((-b0/2),((3*b0)/2),m)
X_axis = ((r_left+r_right)/(b0))

U_right = np.zeros(m)
U_left = np.zeros(m)

for i in range(m):  
    U_right[i] = (circulation_zero * r_right[i]) / ((2 * np.pi) * ((r_right[i])**2 + rc**2))
    U_left[i] = - (circulation_zero * r_left[i]) / ((2 * np.pi) * ((r_left[i])**2 + rc**2))

U_total = (U_left+U_right)/ u_inf

########----->plots
### circulationulation
fig1, ax1 = plt.subplots()

ax1.plot(theta,circulation,".-r", label='circulation')
ax1.set_xlabel('$\\theta$')
ax1.set_ylabel('circulation')
ax1.set_title('circulation Leader')
ax1.legend()

### span loading
fig2, ax2 = plt.subplots()
ax2.plot(eta,SPload, ".-r", label='spanloading')
ax2.set_xlabel('$\\xi$')
ax2.set_ylabel('span loading')
ax2.set_title('Span Loading Leader')
ax2.legend()

### alpha effective
fig3, ax3 = plt.subplots()
ax3.plot((eta),alpha_effe, ".-r", label='effective angle of attack')
ax3.set_xlabel('$\\xi$')
ax3.set_ylabel('effective AOA')
ax3.set_title('effective AOA Leader')
ax3.legend()

### chord distribution
fig4, ax4 = plt.subplots()
ax4.plot(eta,CT,".-r", label='Chord distribution')
ax4.set_xlabel('$\\xi$')
ax4.set_ylabel('Chord distribution')
ax4.set_title('chord distribution Leader')
ax4.legend()

# ### effective Angle of attack
# fig5, ax5 = plt.subplots()
# ax5.plot(eta,alpha_effe,".-r",label='effective AOA')
# ax5.set_xlabel('eta')
# ax5.set_ylabel('Effective AOA')
# ax5.set_title('Effective AOA Leader')
# ax5.legend()

## vortex velocity
fig6, ax6 = plt.subplots()
ax6.plot(X_axis,U_total, ".-r", label='vertical velocity')
ax6.set_xlabel('y/b0')
ax6.set_ylabel('vertical Vortex velocity')
ax6.set_title('vertical velocity Leader')
ax6.legend()


# trailer aircraft (untrimmed)
MATAoddeven = np.zeros((n, n))
MATBoddeven = np.zeros((n, 1))
MATBncoefoddeven = np.zeros((n, 1))
MATA_t = np.zeros((n, n))
MATB_t = np.zeros((n, 1))
MATBncoef_t = np.zeros((n, 1))
MATC_t = np.zeros(n)

N_2 = np.linspace(1,n,n)
U_right1 = np.zeros(m)
U_left1 = np.zeros(m)

Z1 = (a0/(4*b))
Zt = (a0/(4*b))

for i in range(n):
    for j in range(n):
        MATAoddeven[i, j] = np.sin((N_2[j]*np.radians(theta[i]))) + (Z1 * N_2[j] * CT[i] * (np.sin(N_2[j]*np.radians(theta[i])))/(np.sin(np.radians(theta[i]))))

for i in range(n):
    MATBoddeven[i] = 2*Z1*CT[i]

MATBncoefoddeven = np.dot(np.linalg.inv(MATAoddeven), MATBoddeven)

######################################
ds = (b0/2) + (b/2)
r_right1 = np.linspace(ds-(b/2),(ds-(b/2)) + b, m)
r_left1 = np.linspace((ds+b0-(b/2)) , (ds-(b/2))+b0 + b , m)
X_axis1 = ((r_left1+r_right1)/(b0))

for i in range(m):  
    U_right1[i] = (circulation_zero * r_right1[i]) / ((2 * np.pi) * ((r_right1[i])**2 + rc**2))
    U_left1[i] = - (circulation_zero * r_left1[i]) / ((2 * np.pi) * ((r_left1[i])**2 + rc**2))

U_total1= (U_left1+U_right1)/ u_inf

induced_alpha_trailer = U_total1

alpha_wing_trailer = np.zeros(n)
for i in range(n):
    alpha_wing_trailer[i] = (alpha_wing) + induced_alpha_trailer[i]


################################

for i in range(n):
    for j in range(n):
        MATA_t[i, j] = np.sin((N_2[j]*np.radians(theta[i]))) + (Zt * N_2[j] * CT[i] * (np.sin(N_2[j]*np.radians(theta[i])))/(np.sin(np.radians(theta[i]))))

for i in range(n):
    MATB_t[i] = 2*Zt*CT[i]* alpha_wing_trailer[i]

MATBncoef_t = np.dot(np.linalg.inv(MATA_t), MATB_t)

for i in range(n):
    MATC_t[i] = MATBncoef_t[i] - (MATBncoefoddeven[i]*alpha_wing)

#circulation
circulation_t = np.zeros(n)
SUM_t = np.zeros(n)

for i in range(n):
    for j in range(n):
       SUM_t[i] = SUM_t[i] + (MATBncoef_t[j]*np.sin(N_2[j]*np.radians(theta[i])))
    circulation_t[i] = u_inf * b * SUM_t[i]

########----->Q2
alpha_wing_new_t = ((2*Cl)/(np.pi*AR) - MATC_t[0])/MATBncoefoddeven[0]
print("alpha_wing trailer: ",np.rad2deg(alpha_wing_new_t))
######Span Loading Q3

SPload_t = np.zeros(n)
SUM1_t = np.zeros(n)
for i in range(n):
    for j in range(n):
       SUM1_t[i] = SUM1_t[i] + MATBncoef_t[j]*np.sin((N_2[j]*np.radians(theta[i])))
    SPload_t[i] = 2 * AR * SUM1_t[i]

########----->Q4
####calculation of Downwash_angle angle of attack
Downwash_angle_t = np.zeros(n)
SUM2_t = np.zeros(n)
for i in range(n):
    for j in range(n):
            Downwash_angle_t[i] = Downwash_angle_t[i] + (0.5 * (N_2[j] * MATBncoef_t[j] * np.sin(N_2[j] * (np.radians(theta[i]))) / (np.sin((np.radians(theta[i]))))))

####### calculation of alpha effective
alpha_effe_tra = np.zeros(n)
SUM3 = np.zeros(n)
for i in range(n):
    alpha_effe_tra[i] = ((alpha_wing_new_t) + (induced_alpha_trailer[i])) - Downwash_angle_t[i]

#########
addition1 = 0
for i in range(n):
    addition1 = addition1 + circulation_t[i] * (Downwash_angle_t[i] - induced_alpha_trailer[i]) * np.sin(np.radians(theta[i]))* d_theta

drag_t = (0.5*rho*u_inf*b) * addition1
coeff_drag_t = drag_t/ (0.5 * rho * u_inf**2 * s) 
print("CD:",coeff_drag_t)
print("drag_t", drag_t)

###################
addition2 = 0
for i in range(n):
    addition2 = addition2 + circulation_t[i] * np.cos(np.radians(theta[i])) * np.sin(np.radians(theta[i]))* d_theta
moment = (0.25*rho*u_inf*b**2) * addition2
print("Moment: ",moment)
coeff_moment = moment/ (0.5 * rho * u_inf**2 * s * b)
print("CM: ", coeff_moment)

########----->plots

### circulationulation
fig7, ax7 = plt.subplots()

ax7.plot(y/(b/2),circulation,".-r", label='circulation of leader')
ax7.plot((y/(b/2)),circulation_t,".-k", label='circulation of trailer')
ax7.set_xlabel('theta')
ax7.set_ylabel('circulation')
ax7.set_title('circulation Leader vs trailer')
ax7.legend()

### span loading
fig8, ax8 = plt.subplots()
ax8.plot(y/(b/2),SPload, ".-r", label='spanloading of leader')
ax8.plot(y/(b/2),SPload_t, ".-k", label='spanloading of trailer')
ax8.set_xlabel('$\\xi$')
ax8.set_ylabel('span loading')
ax8.set_title('Span Loading Leader vs trailer')
ax8.legend()

### alpha effective
fig9, ax9 = plt.subplots()
ax9.plot((y/(b/2)),alpha_effe, ".-r", label='effective angle of attack of leader')
ax9.plot((y/(b/2)),alpha_effe_tra, ".-k", label='effective angle of attack of trailer')
ax9.set_xlabel('$\\xi$')
ax9.set_ylabel('effective AOA')
ax9.set_title('effective AOA Leader vs trailer')
ax9.legend()

## Induced Angle of attack
fig10, ax10 = plt.subplots()
ax10.plot(X_axis1,induced_alpha_trailer, ".-r", label='induced angle of attack')
ax10.set_xlabel('y/b0')
ax10.set_ylabel('Induced angle of attack')
ax10.set_title('Induced AOA on trailer by Leader')
ax10.legend()

# trailer aircraft (Trimmed)

addition_SUM = 0
summation = 0
D_theo = 0
t = 0.05
f = np.zeros(n)

for i in range(n):
    summation = summation + ((circulation_t[i]) * np.sin(np.radians(theta[i])) * np.cos(np.radians(theta[i])) * (
            np.radians(theta[2]) - np.radians(theta[1])))
addition_SUM = summation

D_theo = (moment * 4) / (rho * u_inf * (b ** 2) * addition_SUM)
print("delta:", np.rad2deg(D_theo))

for i in range(n):
    if 0.7 <= eta[i] <= 0.9:
        f[i] = -1
    elif -0.9 <= eta[i] <= -0.7:
        f[i] = 1
    else:
        f[i] = 0

D_a = D_theo
alpha_t_trim = np.zeros(n)

for i in range(n):
    alpha_t_trim[i] = alpha_wing_new_t + induced_alpha_trailer[i] + (t * D_a * f[i])

MATA3 = np.zeros((n, n))
MATB3 = np.zeros((n, 1))
MATX3 = np.zeros((n, 1))

Z3 = (a0 / (4 * b))

for i in range(n):
    for j in range(n):
        MATA3[i, j] = np.sin((N_2[j] * np.radians(theta[i]))) + (
                (Z3 * N_2[j] * CT[i] * (np.sin(N_2[j] * np.radians(theta[i])))) / (
            np.sin(np.radians(theta[i]))))

for i in range(n):
    MATB3[i] = 2 * Z3 * CT[i] * alpha_t_trim[i]

MATinvA3 = np.linalg.inv(MATA3)
MATX3 = np.dot(MATinvA3, MATB3)

MATD = np.zeros(n)

for i in range(n):
    MATD[i] = MATX3[i] - (MATBncoefoddeven[i] * alpha_wing_new_t) - MATC_t[i]

circulation_trimmed = np.zeros(n)
SUM_TRIM = np.zeros(n)

for i in range(n):
    for j in range(n):
        SUM_TRIM[i] = SUM_TRIM[i] + (MATX3[j]) * np.sin(N_2[j] * np.radians(theta[i]))
    circulation_trimmed[i] = u_inf * b * SUM_TRIM[i]

spanload_trimmed = np.zeros(n)
SUM1_t = np.zeros(n)
for i in range(n):
    for j in range(n):
        SUM1_t[i] = SUM1_t[i] + MATX3[j] * np.sin((N_2[j] * np.radians(theta[i])))
    spanload_trimmed[i] = 2 * AR * SUM1_t[i]

downwash_TRIM = np.zeros(n)
SUM2_TRIM1 = np.zeros(n)
for i in range(n):
    for j in range(n):
        SUM2_TRIM1[i] = SUM2_TRIM1[i] + (N_2[j] * (
                    (MATBncoefoddeven[j] * alpha_wing_new_t) + (MATD[j] * t * D_a) + MATC_t[j]) * (
                                             np.sin(N_2[j] * np.radians(theta[i]))) / (np.sin(np.radians(theta[i]))))
    downwash_TRIM[i] = 0.5 * SUM2_TRIM1[i]

alpha_e_trimmed = np.zeros(n)
for i in range(n):
    alpha_e_trimmed[i] = (alpha_wing_new_t + induced_alpha_trailer[i]) - downwash_TRIM[i]

drag_I_trimmed = np.zeros(n)
integral_TRIM = 0
for i in range(n):
    integral_TRIM = integral_TRIM + circulation_trimmed[i] * (downwash_TRIM[i] - induced_alpha_trailer[i]) * np.sin(
        np.radians(theta[i])) * (np.radians(theta[2]) - np.radians(theta[1]))

drag_I_trimmed = ((rho * u_inf * b) / 2) * integral_TRIM
coefficient_drag_I_trimmed = drag_I_trimmed / (0.5 * rho * u_inf ** 2 * s)

print("Induced drag of trailer: ", drag_I_trimmed)
print("Induced drag coefficient of trailer: ", coefficient_drag_I_trimmed)

########----->plots

### circulationulation
fig11, ax11 = plt.subplots()

ax11.plot(y/(b/2),circulation,".-r", label='circulation of leader')
ax11.plot((y/(b/2)),circulation_trimmed,".-k", label='circulation of trailer')
ax11.set_xlabel('theta')
ax11.set_ylabel('circulation')
ax11.set_title('circulation trimmed case')
ax11.legend()

### span loading
fig12, ax12 = plt.subplots()
ax12.plot(y/(b/2),SPload, ".-r", label='spanloading of leader')
ax12.plot(y/(b/2),spanload_trimmed, ".-k", label='spanloading of trailer')
ax12.set_xlabel('$\\xi$')
ax12.set_ylabel('span loading')
ax12.set_title('Span Loading trimmed case')
ax12.legend()

### alpha effective
fig13, ax13= plt.subplots()
ax13.plot((y/(b/2)),alpha_effe, ".-r", label='effective angle of attack of leader')
ax13.plot((y/(b/2)),alpha_e_trimmed, ".-k", label='effective angle of attack of trailer')
ax13.set_xlabel('$\\xi$')
ax13.set_ylabel('effective AOA')
ax13.set_title('effective AOA trimmed case')
ax13.legend()

### display the figures
plt.show()

