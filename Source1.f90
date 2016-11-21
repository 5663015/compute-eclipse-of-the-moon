      PROGRAM MOON
          implicit none
          REAL *8,DIMENSION(6) :: PO                !输出坐标
          INTEGER :: B                              !目标天体编号
          INTEGER :: C                              !中心天体编号    
          INTEGER :: YEAR=2017                      !年
          INTEGER :: MONTH=8                        !月
          INTEGER :: DAY=7                          !日
          REAL*8 :: DDAY                            !日小数部分
          INTEGER :: HOUR                           !时
          INTEGER :: MIN                            !分
          REAL*8 :: SECOND                          !秒
          REAL *8 :: DJM0                           !儒略日1
          REAL *8 :: DJM                            !儒略日2
          INTEGER :: J                              !状态
          INTEGER :: NDP=3                          !小数秒位数
          INTEGER ,DIMENSION(4) :: IHMSF            !时分秒小数秒数组
          CHARACTER :: SIGN                         !正负号
          CHARACTER*6 :: nams                       !常数名字
          REAL *8::vals                             !常数值
          REAL *8::iteration_val=1.0/(24*3600)      !循环迭代量，1s为多少天
          REAL *8::ERR                              !角度差，判断是否发生食既
          REAL *8::ASUN                             !太阳半径
          REAL *8::RE                               !地球赤道半径
          REAL *8::AM                               !月球半径
          REAL *8::CLIGHT                           !光速
          REAL *8::AU                               !天文单位          
          REAL *8,DIMENSION(3)::E2S                 !地日距离
          REAL *8,DIMENSION(3)::E2O                 !地心与影锥距离
          REAL *8,DIMENSION(3)::E2M                 !地月距离
          REAL *8,DIMENSION(3)::M2O                 !月心与影锥的距离
          REAL *8,DIMENSION(3)::S2M                 !日月矢量
          REAL *8::ANGLE_EOM                        !地锥月角度
          REAL *8::ANGLE_EO                         !地锥角
          REAL *8::ANGLE_MO                         !月锥角
          REAL *8::PRODUCT                          !内积
!************************读出常量**************************
          nams='AU'         !天文单位
          CALL selconQ(nams,vals)
          AU=VALS          
          nams='ASUN'       !太阳半径
          CALL selconQ(nams,vals)
          ASUN=VALS/AU          
          nams='RE'         !地球半径
          CALL selconQ(nams,vals)
          RE=(VALS+65D0)/AU        
          nams='AM'         !月球半径
          CALL selconQ(nams,vals)
          AM=VALS/AU          
          nams='CLIGHT'     !光速
          CALL selconQ(nams,vals)
          CLIGHT=VALS/AU
          
          CALL iau_CAL2JD ( YEAR, MONTH, DAY, DJM0, DJM, J )    !格里历转为儒略日
!***********************循环计算**************************
     DO
         B=11    !太阳
         C=3     !地球
         CALL PLEPH(DJM0+DJM, B, C, PO)       !地日矢量
         CALL PLEPH(DJM0+DJM-((SQRT(PO(1)*PO(1)+PO(2)*PO(2)+PO(3)*PO(3))-ASUN-RE)*iteration_val)/CLIGHT, B, C, PO)!大约8分钟前日地矢量坐标
         E2S(1)=PO(1)                      !地日矢量
         E2S(2)=PO(2)
         E2S(3)=PO(3)         
         E2O(1)=(RE/(RE-ASUN))*E2S(1)      !地锥矢量
         E2O(2)=(RE/(RE-ASUN))*E2S(2)
         E2O(3)=(RE/(RE-ASUN))*E2S(3)

         B=10   !月球
         C=3    !地球
         CALL PLEPH(DJM0+DJM, B, C, PO)
         E2M(1)=PO(1)            !地球到月球矢量
         E2M(2)=PO(2)
         E2M(3)=PO(3)
         M2O(1)=E2O(1)-E2M(1)    !月球到锥点矢量
         M2O(2)=E2O(2)-E2M(2)
         M2O(3)=E2O(3)-E2M(3)
         CALL iau_PDP (E2O,M2O,PRODUCT)      !求内积
         ANGLE_EOM=ACOS(PRODUCT/(SQRT(E2O(1)*E2O(1)+E2O(2)*E2O(2)+E2O(3)*E2O(3))*SQRT(M2O(1)*M2O(1)+M2O(2)*M2O(2)+M2O(3)*M2O(3))))
         ANGLE_EO=ASIN(RE/SQRT(E2O(1)*E2O(1)+E2O(2)*E2O(2)+E2O(3)*E2O(3)))
         ANGLE_MO=ASIN(AM/SQRT(M2O(1)*M2O(1)+M2O(2)*M2O(2)+M2O(3)*M2O(3)))
         ERR=ANGLE_EOM-ANGLE_EO-ANGLE_MO       !求角度之差
         DJM=DJM+iteration_val    !儒略日增加1s
         IF(ERR.LE.0D0) EXIT      !此刻ERR小于0且上时刻ERR_LASE大于0时退出         
     END DO
!*************************结果输出******************************     
         CALL iau_JD2CAL ( DJM0, DJM, YEAR, MONTH, DAY, DDAY, J )     !儒略日转为格里历
         DDAY=DDAY-(32.164+37)*iteration_val               !转为UCT
         CALL sla_DD2TF (NDP, DDAY, SIGN, IHMSF)           !小数天转为时分秒小数秒.
         WRITE (*,*)'初亏时间：'
         WRITE (*,*) YEAR,'年',MONTH,'月',DAY,'日',IHMSF(1),':',IHMSF(2),':',IHMSF(3)+IHMSF(4)/1000.0   !输出结果
        
      END PROGRAM MOON

