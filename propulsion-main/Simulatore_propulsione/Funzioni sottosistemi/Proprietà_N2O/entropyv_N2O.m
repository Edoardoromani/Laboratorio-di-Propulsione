function[sv]=entropyv_N2O(Told)
%Temperature (K)
T=[183.00
184.00
185.00
186.00
187.00
188.00
189.00
190.00
191.00
192.00
193.00
194.00
195.00
196.00
197.00
198.00
199.00
200.00
201.00
202.00
203.00
204.00
205.00
206.00
207.00
208.00
209.00
210.00
211.00
212.00
213.00
214.00
215.00
216.00
217.00
218.00
219.00
220.00
221.00
222.00
223.00
224.00
225.00
226.00
227.00
228.00
229.00
230.00
231.00
232.00
233.00
234.00
235.00
236.00
237.00
238.00
239.00
240.00
241.00
242.00
243.00
244.00
245.00
246.00
247.00
248.00
249.00
250.00
251.00
252.00
253.00
254.00
255.00
256.00
257.00
258.00
259.00
260.00
261.00
262.00
263.00
264.00
265.00
266.00
267.00
268.00
269.00
270.00
271.00
272.00
273.00
274.00
275.00
276.00
277.00
278.00
279.00
280.00
281.00
282.00
283.00
284.00
285.00
286.00
287.00
288.00
289.00
290.00
291.00
292.00
293.00
294.00
295.00
296.00
297.00
298.00
299.00
300.00
301.00
302.00
303.00
304.00
305.00
306.00
307.00
308.00
309.00];

%Vapor Entropy (kJ/kg-K)
s=[2.0399
    2.0320
    2.0242
    2.0165
    2.0089
    2.0015
    1.9941
    1.9868
    1.9797
    1.9726
    1.9656
    1.9587
    1.9519
    1.9452
    1.9386
    1.9321
    1.9256
    1.9193
    1.9130
    1.9068
    1.9006
    1.8945
    1.8885
    1.8826
    1.8768
    1.8710
    1.8652
    1.8596
    1.8540
    1.8484
    1.8429
    1.8375
    1.8321
    1.8268
    1.8215
    1.8163
    1.8111
    1.8060
    1.8009
    1.7959
    1.7909
    1.7860
    1.7811
    1.7762
    1.7714
    1.7666
    1.7618
    1.7571
    1.7524
    1.7478
    1.7431
    1.7385
    1.7340
    1.7294
    1.7249
    1.7204
    1.7160
    1.7115
    1.7071
    1.7027
    1.6983
    1.6939
    1.6896
    1.6852
    1.6809
    1.6766
    1.6722
    1.6679
    1.6636
    1.6594
    1.6551
    1.6508
    1.6465
    1.6422
    1.6379
    1.6336
    1.6293
    1.6250
    1.6207
    1.6164
    1.6121
    1.6077
    1.6033
    1.5989
    1.5945
    1.5901
    1.5856
    1.5811
    1.5766
    1.5721
    1.5675
    1.5628
    1.5581
    1.5534
    1.5486
    1.5438
    1.5389
    1.5339
    1.5289
    1.5237
    1.5185
    1.5132
    1.5078
    1.5023
    1.4967
    1.4910
    1.4851
    1.4790
    1.4728
    1.4664
    1.4598
    1.4530
    1.4459
    1.4385
    1.4309
    1.4228
    1.4143
    1.4053
    1.3957
    1.3854
    1.3742
    1.3619
    1.3480
    1.3321
    1.3129
    1.2877
    1.2466];
sv=interp1(T,s,Told)*1000;