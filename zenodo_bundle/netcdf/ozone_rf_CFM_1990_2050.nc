CDF      
      time       lats            CDI       ?Climate Data Interface version 1.6.8 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.4     history      FThu Dec 17 11:39:13 2020: /home/opt/package/nco/linux64/bin/ncrcat ozone_rf_1990_2010_sngl.nc ozone_rf_SLCF_2015_2050.nc ozone_rf_SLCF_1990_2050.nc
Thu Dec 17 11:34:51 2020: /home/opt/package/nco/linux64/bin/ncks -v time,lats,drfrc,dteqdr tst5.nc tst6.nc
Thu Dec 17 11:34:44 2020: /home/opt/package/nco/linux64/bin/ncap2 -s dteqdr=dteqdrt tst4.nc tst5.nc
Thu Dec 17 11:34:31 2020: /home/opt/package/nco/linux64/bin/ncap2 -s drfrc=drfrct tst3.nc tst4.nc
Thu Dec 17 11:34:16 2020: /home/opt/package/nco/linux64/bin/ncks -v time,lats,drfrct,dteqdrt tst2.nc tst3.nc
Thu Dec 17 11:34:10 2020: /home/opt/package/nco/linux64/bin/ncap2 -s drfrct=drfrc.total($region) tst1.nc tst2.nc
Thu Dec 17 11:33:57 2020: /home/opt/package/nco/linux64/bin/ncap2 -s dteqdrt=dteqdr.total($region) tst.nc tst1.nc
Thu Dec 17 11:33:36 2020: /home/opt/package/nco/linux64/bin/ncks -d region,0 ozone_rf_1990_2010_corr.nc tst.nc
Tue Aug 18 13:43:00 2020: cdo mulc,-1 ozone_rf_1990_2010.nc ozone_rf_1990_2010_corr.nc
Fri Aug 14 10:00:22 2020: cdo selyear,1990,1995,2000,2005,2010 ozone_rf_1990_2015.nc ozone_rf_1990_2010.nc      CDO       CClimate Data Operators version 1.6.8rc2 (http://mpimet.mpg.de/cdo)     NCO       4.4.7      nco_openmp_thread_number                  drfrc                      
_FillValue        �Ç�       	long_name         TOA direct radiative forcing   missing_value         �Ç�       units         W m-2            �   dteqdr                     
_FillValue        �Ç�       	long_name         <Equilibrium temperature change from direct radiative forcing   missing_value         �Ç�       units         K            �   lats               	long_name         Receptor region    units         .1 - Arctic, 2 - NH midlat, 3 - Tropics, 4 - SH     axis      Y            �   time                standard_name         time   	long_name         years      units         year as %Y.%f      calendar      proleptic_gregorian         �?�      @       @      @      ?��B_�'I?��6���J���#s6���8�gH�}c*:�>��|�"��n���˘%ƿ�j���@�     ?���9�bS?�)=w������Ĳ����Fxu,ۿmc���VԿUU���2����S�P����=p�@�,     ?�kua�t�?��2�ꅿ��{���a���˦\r׿n�.^��e?��
�N��;�˶��G/l?@�@     ?�� ��?�[�cAY�����ʿ�5�����f�=G�>��_��R,�l�{�G��|�$W�'@�T     ?tG����_?a�3k�	��w�z=MN��t���S�a�I�
�b���}�D�s�f����p5;99�@�h                                                                     @�|     �v�B����0'#��c?lM�a�(?q$��u~���\���d�U�p��L?8�F��{?_�?ˁ@��     ��?�sy ۿ��A�_V3���$ʒn�����DǼ��W/^�%��`�rဎ��.qf���kr����@��     ���}��e����%o���Qm?��ĕ��W����,ץJ���K��T����k��yT�\�� @��     ���F�v����`�~����
�?�����Ó1������翥j�L7���5(]���Bۂ*@��     ������&��#�ZI����:�;��������@����N�#��JO(۳��=g]׷���S��+@�     