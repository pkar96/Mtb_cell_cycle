classdef Cell < matlab.mixin.Copyable
   properties
      t
      v
      rate
      oris
      vNextDivs = [];
      tNextDivs = [];
      vNextInit
      vOfInits
      tOfInits
      oOfInits
      vb
      vi
      vd
      tLastDiv
      tLastInit
      tNextDiv_in
      vNextDiv_in
      vNextDiv_dv
      tNextDiv_dv
      divr
      moth_info
   end
   methods
      function obj = Cell(v)
         if nargin == 0
            obj.t = 0;
            obj.v = 1;
            obj.vb = 1;
            obj.vi = 1;
            obj.vd = 1;
            obj.rate = 0;
            obj.vNextInit = [];
            obj.oris = [];
            obj.vNextDivs = [];
            obj.tNextDivs = [];
            obj.tLastDiv = 0;
            obj.tLastInit = 0;
            obj.vOfInits = [1];
            obj.tOfInits = [0];
            obj.oOfInits = [1];
            obj.tNextDiv_in=[];
            obj.vNextDiv_in=[];
            obj.vNextDiv_dv=NaN;
            obj.tNextDiv_dv=NaN;
            obj.divr = 0.5;
            obj.moth_info=NaN;
         else
            obj.t = 0;
            obj.v = v;
            obj.vb = v;
            obj.vi = v;
            obj.vd = v;
            obj.rate = 0;
            obj.vNextInit = [];
            obj.oris = [];
            obj.vNextDivs = [];
            obj.tNextDivs = [];
            obj.tLastDiv = 0;
            obj.tLastInit = 0;
            obj.vOfInits = [v];
            obj.tOfInits = [0];
            obj.oOfInits = [2];
            obj.tNextDiv_in=[];
            obj.vNextDiv_in=[];
            obj.vNextDiv_dv=NaN;
            obj.tNextDiv_dv=NaN;
            obj.divr = 0.5;
            obj.moth_info=NaN;
         end
      end
   end
end