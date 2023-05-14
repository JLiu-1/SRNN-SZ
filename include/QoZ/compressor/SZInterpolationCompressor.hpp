#ifndef _SZ_INTERPOLATION_COMPRESSSOR_HPP
#define _SZ_INTERPOLATION_COMPRESSSOR_HPP
#include "QoZ/compressor/Compressor.hpp"
#include "QoZ/predictor/Predictor.hpp"
#include "QoZ/predictor/LorenzoPredictor.hpp"
#include "QoZ/quantizer/Quantizer.hpp"
#include "QoZ/encoder/Encoder.hpp"
#include "QoZ/lossless/Lossless.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/utils/MemoryUtil.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/utils/FileUtil.hpp"
#include "QoZ/utils/Interpolators.hpp"
#include "QoZ/utils/Timer.hpp"
#include "QoZ/def.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/utils/CoeffRegression.hpp"
#include "QoZ/utils/Sample.hpp"
#include <cstring>
#include <cmath>
#include <limits>
namespace QoZ {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZInterpolationCompressor : public concepts::CompressorInterface<T> {//added heritage
    public:


        SZInterpolationCompressor(Quantizer quantizer, Encoder encoder, Lossless lossless) :
                quantizer(quantizer), encoder(encoder), lossless(lossless) {

            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) {
            size_t remaining_length = cmpSize;
            uchar *buffer = lossless.decompress(cmpData, remaining_length);
            int levelwise_predictor_levels;
            bool blockwiseTuning;
            uchar const *buffer_pos = buffer;


            


            std::vector <uint8_t> interpAlgo_list;
            std::vector <uint8_t> interpDirection_list;
            std::vector <uint8_t> cubicSplineType_list;
            


            std::vector <QoZ::Interp_Meta> interpMeta_list;
            int fixBlockSize;
            int trimToZero;
           
            read(global_dimensions.data(), N, buffer_pos, remaining_length);        
            read(blocksize, buffer_pos, remaining_length);
            /*
            read(interpolator_id, buffer_pos, remaining_length);           
            read(direction_sequence_id, buffer_pos, remaining_length);  
            read(cubicSplineType, buffer_pos, remaining_length);         
            */
            read(interp_meta, buffer_pos, remaining_length);    
            read(alpha,buffer_pos,remaining_length);
            read(beta,buffer_pos,remaining_length);
            read(maxStep,buffer_pos,remaining_length);
           
            read(levelwise_predictor_levels,buffer_pos, remaining_length);
            read(blockwiseTuning,buffer_pos, remaining_length);
            read(fixBlockSize,buffer_pos, remaining_length);
            int fused_dim=-1;
            read(fused_dim,buffer_pos, remaining_length);
            
           // size_t cross_block=0;
            //read(cross_block,buffer_pos, remaining_length);
            read(trimToZero,buffer_pos, remaining_length);
            //int blockOrder=0;
            //read(blockOrder,buffer_pos, remaining_length); 
            //int regressiveInterp;   
            //read(regressiveInterp,buffer_pos, remaining_length);       
            if (trimToZero>0){
                quantizer.setTrimToZero(trimToZero);
            }
           
            if(blockwiseTuning){
                size_t meta_num;
                read(meta_num,buffer_pos, remaining_length);
                interpMeta_list.resize(meta_num);
                read(interpMeta_list.data(),meta_num,buffer_pos, remaining_length);
            }

            else if(levelwise_predictor_levels>0){
                /*
                interpAlgo_list=std::vector <uint8_t>(levelwise_predictor_levels,0);
                interpDirection_list=std::vector <uint8_t>(levelwise_predictor_levels,0);
                cubicSplineType_list.resize(levelwise_predictor_levels);
                read(interpAlgo_list.data(),levelwise_predictor_levels,buffer_pos, remaining_length);
                read(interpDirection_list.data(),levelwise_predictor_levels,buffer_pos, remaining_length);
                read(cubicSplineType_list.data(),levelwise_predictor_levels,buffer_pos, remaining_length);
                */
                interpMeta_list.resize(levelwise_predictor_levels);
                read(interpMeta_list.data(),levelwise_predictor_levels,buffer_pos, remaining_length);
                //for(auto meta:interpMeta_list)
                //    QoZ::print_meta(meta);

            }           
            init();   
          
            //QoZ::Timer timer(true);
            quantizer.load(buffer_pos, remaining_length);
            encoder.load(buffer_pos, remaining_length);
            quant_inds = encoder.decode(buffer_pos, num_elements);
            encoder.postprocess_decode();

            lossless.postdecompress_data(buffer);
            //timer.stop("decode");
            //timer.start();
            double eb = quantizer.get_eb();
            if(!anchor){
                *decData = quantizer.recover(0, quant_inds[quant_index++]);
            }
            
            else{
                recover_grid(decData,global_dimensions,maxStep,fused_dim);                   
                interpolation_level--;           
            }
            size_t meta_index=0;
            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {

                if (alpha<0) {
                    if (level >= 3) {
                        quantizer.set_eb(eb * eb_ratio);
                    } else {
                        quantizer.set_eb(eb);
                    }
                }
                
                else if (alpha>=1){
                    
                    
                    double cur_ratio=pow(alpha,level-1);
                    if (cur_ratio>beta){
                        cur_ratio=beta;
                    }
                    
                    quantizer.set_eb(eb/cur_ratio);
                }
                else{
                    
                    
                    double cur_ratio=1-(level-1)*alpha;
                    if (cur_ratio<beta){
                        cur_ratio=beta;
                    }
                   
                    quantizer.set_eb(eb*cur_ratio);
                }
                /*
                uint cur_interpolator=interpolator_id;
                uint cur_direction=direction_sequence_id;
                uint cur_splinetype=cubicSplineType;
                */
                QoZ::Interp_Meta cur_meta;
                if(blockwiseTuning)
                    cur_meta=interpMeta_list[meta_index++];
                else if (levelwise_predictor_levels==0){
                    cur_meta=interp_meta;
                }
                else{

                    if (level-1<levelwise_predictor_levels){
                        /*
                        cur_interpolator=interpAlgo_list[level-1];
                        cur_direction=interpDirection_list[level-1];
                        cur_splinetype=cubicSplineType_list[level-1];
                        */
                        cur_meta=interpMeta_list[level-1];
                    }
                    else{
                        /*
                        cur_interpolator=interpAlgo_list[levelwise_predictor_levels-1];
                        cur_direction=interpDirection_list[levelwise_predictor_levels-1];
                        cur_splinetype=cubicSplineType_list[levelwise_predictor_levels-1];
                        */
                        cur_meta=interpMeta_list[levelwise_predictor_levels-1];
                    }
                }
                     
                size_t stride = 1U << (level - 1);
                size_t cur_blocksize;
                
                if (fixBlockSize>0){
                    cur_blocksize=fixBlockSize;
                }
                else{
                    cur_blocksize=blocksize*stride;
                }
                
                auto inter_block_range = std::make_shared<QoZ::multi_dimensional_range<T, N>>(decData,std::begin(global_dimensions), std::end(global_dimensions),
                                                           cur_blocksize, 0,0);//blockOrder);
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                
                //timer.stop("prep");
                
                for (auto block = inter_begin; block != inter_end; ++block) {

                    auto start_idx=block.get_global_index();
                    auto end_idx = start_idx;


                    for (int i = 0; i < N; i++) {
                        end_idx[i] += cur_blocksize;
                        if (end_idx[i] > global_dimensions[i] - 1) {
                            end_idx[i] = global_dimensions[i] - 1;
                        }
                    }
               
                     block_interpolation(decData, block.get_global_index(), end_idx, PB_recover,
                                        interpolators[cur_meta.interpAlgo], cur_meta, stride,0,0,0);//,cross_block,regressiveInterp);
                
                        


                

                }
               
            }
            quantizer.postdecompress_data();
            return decData;
        }
        
        T *decompress_block(uchar const *cmpData, const size_t &cmpSize, T *decData) {
            /*
            size_t remaining_length = cmpSize;
            uchar *buffer = lossless.decompress(cmpData, remaining_length);
            int levelwise_predictor_levels;
            bool blockwiseTuning;
            uchar const *buffer_pos = buffer;
            std::vector <uint8_t> interpAlgo_list;
            std::vector <uint8_t> interpDirection_list;
            int fixBlockSize;
            int trimToZero;

            read(global_dimensions.data(), N, buffer_pos, remaining_length);    
            read(blocksize, buffer_pos, remaining_length);
            read(interpolator_id, buffer_pos, remaining_length);
            read(direction_sequence_id, buffer_pos, remaining_length);
            read(alpha,buffer_pos,remaining_length);
            read(beta,buffer_pos,remaining_length);
            read(maxStep,buffer_pos,remaining_length);
            read(levelwise_predictor_levels,buffer_pos, remaining_length);
            read(blockwiseTuning,buffer_pos, remaining_length);
            read(fixBlockSize,buffer_pos, remaining_length);

            size_t cross_block=0;
            read(cross_block,buffer_pos, remaining_length);
            read(trimToZero,buffer_pos, remaining_length);
            int blockOrder=0;
            read(blockOrder,buffer_pos, remaining_length);
            if (trimToZero>0){
                quantizer.setTrimToZero(trimToZero);
            }

            if(blockwiseTuning){
                size_t ops_num;
                read(ops_num,buffer_pos, remaining_length);
                interpAlgo_list=std::vector <uint8_t>(ops_num,0);
                interpDirection_list=std::vector <uint8_t>(ops_num,0);
                read(interpAlgo_list.data(),ops_num,buffer_pos, remaining_length);
                read(interpDirection_list.data(),ops_num,buffer_pos, remaining_length);
            }
            
            else if(levelwise_predictor_levels>0){
                interpAlgo_list=std::vector <uint8_t>(levelwise_predictor_levels,0);
                interpDirection_list=std::vector <uint8_t>(levelwise_predictor_levels,0);
                read(interpAlgo_list.data(),levelwise_predictor_levels,buffer_pos, remaining_length);
                read(interpDirection_list.data(),levelwise_predictor_levels,buffer_pos, remaining_length);
            }

            init();
            //QoZ::Timer timer(true);
            quantizer.load(buffer_pos, remaining_length);
            encoder.load(buffer_pos, remaining_length);
            quant_inds = encoder.decode(buffer_pos, num_elements);

            encoder.postprocess_decode();

            lossless.postdecompress_data(buffer);
            //timer.stop("decode");
            //timer.start();
            double eb = quantizer.get_eb();
            if(!anchor){
                *decData = quantizer.recover(0, quant_inds[quant_index++]);
            }
            
            else{
                recover_grid(decData,global_dimensions,maxStep);  
                interpolation_level--;     
            }
            size_t op_index=0;

            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
                if (alpha<0) {
                    if (level >= 3) {
                        quantizer.set_eb(eb * eb_ratio);
                    } else {
                        quantizer.set_eb(eb);
                    }
                }
                else if (alpha>=1){
                    double cur_ratio=pow(alpha,level-1);
                    if (cur_ratio>beta){
                        cur_ratio=beta;
                    }
                    quantizer.set_eb(eb/cur_ratio);
                }
                else{
                    double cur_ratio=1-(level-1)*alpha;
                    if (cur_ratio<beta){
                        cur_ratio=beta;
                    }
                    quantizer.set_eb(eb*cur_ratio);
                }
                uint8_t cur_interpolator=interpolator_id;
                uint8_t cur_direction=direction_sequence_id;
                
                size_t stride = 1U << (level - 1);
                size_t cur_blocksize=blocksize;
            
                auto inter_block_range = std::make_shared<
                        QoZ::multi_dimensional_range<T, N>>(decData,
                                                           std::begin(global_dimensions), std::end(global_dimensions),
                                                           cur_blocksize, 0,blockOrder);
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                
                //timer.stop("prep");
                    for (auto block = inter_begin; block != inter_end; ++block) {
                        auto start_idx=block.get_global_index();
                        auto end_idx = start_idx;
                        for (int i = 0; i < N; i++) {
                            end_idx[i] += cur_blocksize;
                            if (end_idx[i] > global_dimensions[i] - 1) {
                                end_idx[i] = global_dimensions[i] - 1;
                            }
                        }
                        block_interpolation_block3d(decData, start_idx, end_idx, PB_recover,
                                            interpolators[interpAlgo_list[op_index]], interpDirection_list[op_index], stride,0,cross_block);
                        op_index++;
                    }    
            }
            quantizer.postdecompress_data();
            */
            return decData;
        }
        
        uchar *compress(Config &conf, T *data, size_t &compressed_size,int tuning=0) {
            return compress(conf,data,compressed_size,tuning,0,0);
        }
        // compress given the error bound
        uchar *compress( Config &conf, T *data, size_t &compressed_size,int tuning,int start_level,int end_level=0) {
            //tuning 0: normal compress 1:tuning to return qbins and psnr 2: tuning to return prediction loss
            Timer timer;
            timer.start();
            if (conf.trimToZero>0){
                quantizer.setTrimToZero(conf.trimToZero);
            }
            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
            blocksize = conf.interpBlockSize;
            maxStep=conf.maxStep;
            /*
            interpolator_id = conf.interpMeta.interpAlgo;
            interp_paradigm=conf.interpMeta.interpParadigm;
            cubicSplineType=conf.cubicSplineType;
            direction_sequence_id = conf.interpMeta.interpDirection;
            adj_interp = conf.interpMeta.adjInterp;
            */
            interp_meta=conf.interpMeta;

            alpha=conf.alpha;
            beta=conf.beta;
            std::vector<Interp_Meta>interp_metas;
            //size_t cross_block=conf.crossBlock;
            //int regressiveInterp=conf.regressiveInterp;
            init();
            if (tuning){
                std::vector<int>().swap(quant_inds);
                std::vector<int>().swap(conf.quant_bins);
                conf.quant_bin_counts=std::vector<size_t>(interpolation_level,0);
                conf.decomp_square_error=0.0;

            }
            /*
            if(tuning==0 and conf.peTracking){
                prediction_errors.resize(num_elements,0);
                peTracking=1;
            }
            */
            quant_inds.reserve(num_elements);
            size_t interp_compressed_size = 0;
            double eb = quantizer.get_eb();

            if (start_level<=0 or start_level>interpolation_level ){

                start_level=interpolation_level;

                
            } 
            if(end_level>=start_level or end_level<0){
                end_level=0;
            }


            if(!anchor){
                quant_inds.push_back(quantizer.quantize_and_overwrite(*data, 0));
            }
            else if (start_level==interpolation_level){
                if(tuning){
                    conf.quant_bin_counts[start_level-1]=quant_inds.size();
                }
                build_grid(conf,data,maxStep,tuning);
                start_level--;
            }
            double predict_error=0.0;
            int levelwise_predictor_levels=conf.interpMeta_list.size();

            for (uint level = start_level; level > end_level && level <= start_level; level--) {
                cur_level=level;
                if (alpha<0) {
                    if (level >= 3) {
                        quantizer.set_eb(eb * eb_ratio);
                    } else {
                        quantizer.set_eb(eb);
                    }
                }
                else if (alpha>=1){              
                    double cur_ratio=pow(alpha,level-1);
                    if (cur_ratio>beta){
                        cur_ratio=beta;
                    }            
                    quantizer.set_eb(eb/cur_ratio);
                }
                else{              
                    double cur_ratio=1-(level-1)*alpha;
                    if (cur_ratio<beta){
                        cur_ratio=beta;
                    }             
                    quantizer.set_eb(eb*cur_ratio);
                }
                /*
                uint8_t cur_interpolator;
                uint8_t cur_paradigm;
                uint8_t cur_splinetype;
                uint8_t cur_direction;
                uint8_t cur_adj;
                */
                QoZ::Interp_Meta cur_meta;
                if (levelwise_predictor_levels==0){
                    cur_meta=interp_meta;
                    
                }
                else{
                    if (level-1<levelwise_predictor_levels){
                        cur_meta=conf.interpMeta_list[level-1];
                    }
                    else{
                        cur_meta=conf.interpMeta_list[levelwise_predictor_levels-1];
                    }
                    /*
                    if(level==1 and conf.adaptiveMultiDimStride>0 and tuning<=0){
                        std::vector<double> vars;
                        QoZ::calculate_interp_error_vars<T,N>(data,  conf.dims,vars,cur_meta.interpAlgo,cur_meta.cubicSplineType,conf.adaptiveMultiDimStride,0);
                        QoZ::preprocess_vars<N>(vars);
                        for(size_t i=0;i<N;i++)
                            cur_meta.dimCoeffs[i]=vars[i];
                        conf.interpMeta_list[0]=cur_meta;
                    }
                    */
                }
                
                size_t stride = 1U << (level - 1);
                size_t cur_blocksize;
                if (conf.fixBlockSize>0){
                    cur_blocksize=conf.fixBlockSize;
                }
                else{
                    cur_blocksize=blocksize*stride;
                }       
                auto inter_block_range = std::make_shared<
                        QoZ::multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                           std::end(global_dimensions),
                                                           cur_blocksize, 0,0);//conf.blockOrder);
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {
                    auto start_idx=block.get_global_index();
                    auto end_idx = start_idx;
                    for (int i = 0; i < N; i++) {
                        end_idx[i] += cur_blocksize ;
                        if (end_idx[i] > global_dimensions[i] - 1) {
                            end_idx[i] = global_dimensions[i] - 1;
                        }
                    }
                    if(!conf.blockwiseTuning){
                        /*
                        if(peTracking)
                            
                            predict_error+=block_interpolation(data, start_idx, end_idx, PB_predict_overwrite,
                                        interpolators[cur_meta.interpAlgo],cur_meta, stride,3,0,0);//,cross_block,regressiveInterp);

                        else */
                            predict_error+=block_interpolation(data, start_idx, end_idx, PB_predict_overwrite,
                                        interpolators[cur_meta.interpAlgo],cur_meta, stride,tuning,0,0);//,cross_block,regressiveInterp);

                    }
                    else{
                        size_t min_len=16;
                        auto start_idx=block.get_global_index();
                        auto end_idx = start_idx;
                        //std::array<size_t,N> block_lengths;
                        std::array<size_t,N> sample_starts;
                        //std::array<size_t,N> sample_strides;
                        for (int i = 0; i < N; i++) {
                            end_idx[i] += cur_blocksize ;
                            if (end_idx[i] > global_dimensions[i] - 1) {
                                end_idx[i] = global_dimensions[i] - 1;
                            }
                            double cur_rate=conf.blockwiseSampleRate;
                            size_t  cur_length=(end_idx[i]-start_idx[i])+1,cur_stride=stride*cur_rate;
                            while(cur_stride>stride){
                                if(cur_length/cur_stride>=min_len)
                                    break;
                                cur_stride/=2;
                                cur_rate/=2;
                                if(cur_stride<stride){
                                    cur_stride=stride;
                                    cur_rate=1;
                                }
                            }
                            //sample_strides[i]=cur_stride;
                            double temp=0.5-0.5/cur_rate;
                            sample_starts[i]=temp*cur_length+start_idx[i];
                            std::cout<<start_idx[i]<<" "<<end_idx[i]<<" "<<sample_starts[i]<<" "<<stride<<std::endl;

                        }
                         std::cout<<"----"<<std::endl;
                        std::vector<T> orig_sampled_block;
                        size_t local_idx=0;
                        std::array<size_t,N> sb_starts;
                        std::fill(sb_starts.begin(),sb_starts.end(),0);
                        std::array<size_t,N> sb_dims;
                        std::fill(sb_dims.begin(),sb_dims.end(),0);
                        size_t x,y,z;
                        if(N==2){
                            for(size_t x=sample_starts[0];x<=end_idx[0] ;x+=stride){
                                sb_dims[0]++;
                                for(size_t y=sample_starts[1];y<=end_idx[1];y+=stride){
                                    sb_dims[1]++;
                                    size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1];
                                    orig_sampled_block.push_back(data[global_idx]);
                                    
                                }
                            }
                        }
                        else if(N==3){
                            for(size_t x=sample_starts[0];x<=end_idx[0] ;x+=stride){
                                sb_dims[0]++;
                                for(size_t y=sample_starts[1];y<=end_idx[1];y+=stride){
                                    sb_dims[1]++;
                                    for(size_t z=sample_starts[2];z<=end_idx[2];z+=stride){
                                        sb_dims[2]++;
                                        size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1]+z*dimension_offsets[2];
                                        orig_sampled_block.push_back(data[global_idx]);
                                    }
                                }
                            }
                        } 

                        QoZ::Interp_Meta best_meta,cur_meta;
                        double best_loss=std::numeric_limits<double>::max();
                        std::vector<uint8_t> interpAlgo_Candidates={QoZ::INTERP_ALGO_LINEAR, QoZ::INTERP_ALGO_CUBIC};
                        std::vector<uint8_t> interpParadigm_Candidates={0};
                        std::vector<uint8_t> cubicSplineType_Candidates={0};
                        std::vector<uint8_t> interpDirection_Candidates={0, QoZ::factorial(N) -1};
                        if(conf.fused_dim>=0){
                            if(conf.fused_dim==0)
                                interpDirection_Candidates={6,7};
                            else if (conf.fused_dim==1)
                                interpDirection_Candidates={8,9};
                            else
                                interpDirection_Candidates={10,11};
                        }
                        std::vector<uint8_t> adjInterp_Candidates={0};

                        if(conf.multiDimInterp>0){
                            for(size_t i=1;i<=conf.multiDimInterp;i++)
                                interpParadigm_Candidates.push_back(i);
                       
                        }

                        if (conf.naturalSpline){
                            cubicSplineType_Candidates.push_back(1);
                        }
                        if(conf.fullAdjacentInterp){
                            adjInterp_Candidates.push_back(1);
                            //for(size_t i=1;i<=conf.fullAdjacentInterp;i++)
                            //    adjInterp_Candidates.push_back(i);
                        }

                        std::vector<T> cur_block;
                        for (auto &interp_op: interpAlgo_Candidates) {
                            cur_meta.interpAlgo=interp_op;
                            for (auto &interp_pd: interpParadigm_Candidates) {
                                if(conf.fused_dim>=0 and interp_pd>1)
                                    continue;
                                cur_meta.interpParadigm=interp_pd;

                                for (auto &interp_direction: interpDirection_Candidates) {
                                    if (conf.fused_dim<0 and (interp_pd==1 or  (interp_pd==2 and N<=2)) and interp_direction!=0)
                                        continue;
                                    cur_meta.interpDirection=interp_direction;
                                    for(auto &cubic_spline_type:cubicSplineType_Candidates){
                                        if (interp_op!=QoZ::INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                            break;
                                        cur_meta.cubicSplineType=cubic_spline_type;
                                        for(auto adj_interp:adjInterp_Candidates){
                                            if (interp_op!=QoZ::INTERP_ALGO_CUBIC and adj_interp!=0)
                                                break;
                                        
                                            cur_meta.adjInterp=adj_interp;
                                            cur_block=orig_sampled_block;
                                            
                                            double cur_loss=block_interpolation(cur_block.data(), sb_starts, sb_dims, PB_predict_overwrite,
                                                interpolators[cur_meta.interpAlgo],cur_meta, 1,2,0,0);//,cross_block,regressiveInterp);

                                            if(cur_loss<best_loss){
                                                best_loss=cur_loss;
                                                best_meta=cur_meta;
                                            }

                                            size_t local_idx=0;
                                        }
                                    }
                                }
                            }
                        }
                        interp_metas.push_back(best_meta);
                        predict_error+=block_interpolation(data, start_idx, end_idx, PB_predict_overwrite,
                                        interpolators[best_meta.interpAlgo],best_meta, stride,tuning,0,0);//,cross_block,regressiveInterp);
                    }

                    
                        
                }
                if(tuning){
                    conf.quant_bin_counts[level-1]=quant_inds.size();
                }
            }                    
            //timer.start();

            quantizer.set_eb(eb);
            /*
            if(peTracking){
                conf.predictionErrors=prediction_errors;
            }
            */
            if (tuning){
                conf.quant_bins=quant_inds;
                std::vector<int>().swap(quant_inds);
                conf.decomp_square_error=predict_error;
                size_t bufferSize = 1;
                uchar *buffer = new uchar[bufferSize];
                buffer[0]=0;
                return buffer;
            }
            /*
            if(peTracking){
                QoZ::writefile<float>("interp_pred.errors", prediction_errors.data(),prediction_errors.size());//added.

            }
            */
            if(conf.verbose)
                timer.stop("prediction");
            /*
            for(size_t i=0;i<num_elements;i++){
                if(!mark[i])
                    std::cout<<i<<std::endl;
            }
             */
            
            
            //timer.start();
            assert(quant_inds.size() == num_elements);
            encoder.preprocess_encode(quant_inds, 0);
            size_t bufferSize = 1.2 * (quantizer.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size());
            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;
            write(global_dimensions.data(), N, buffer_pos);
            write(blocksize, buffer_pos);
            /*
            write(interp_meta.interpAlgo, buffer_pos);
            write(interp_meta.interpParadigm, buffer_pos);
            write(interp_meta.cubicSplineType, buffer_pos);
            write(interp_meta.interpDirection, buffer_pos);
            write(interp_meta.adjInterp, buffer_pos);
            */
            write(interp_meta, buffer_pos);
            
            write(alpha,buffer_pos);
            write(beta,buffer_pos);
            write(maxStep,buffer_pos);
            write(levelwise_predictor_levels,buffer_pos);
            write(conf.blockwiseTuning,buffer_pos);
            write(conf.fixBlockSize,buffer_pos);
            write(conf.fused_dim,buffer_pos);
            //write(cross_block,buffer_pos);
            write(conf.trimToZero,buffer_pos);
            //write(conf.blockOrder,buffer_pos);
            //write(conf.regressiveInterp,buffer_pos);
            if(conf.blockwiseTuning){
                size_t meta_num=interp_metas.size();
                write(meta_num,buffer_pos);
                write(interp_metas.data(),meta_num,buffer_pos);

            }
            else if(levelwise_predictor_levels>0){
                write(conf.interpMeta_list.data(),levelwise_predictor_levels,buffer_pos);
               
            }
            quantizer.save(buffer_pos);
            quantizer.postcompress_data();
            quantizer.clear();
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();          
            //timer.stop("Coding");
            //timer.start();
            assert(buffer_pos - buffer < bufferSize);         
            uchar *lossless_data = lossless.compress(buffer,
                                                     buffer_pos - buffer,
                                                     compressed_size);
            lossless.postcompress_data(buffer);
            //timer.stop("Lossless") ;
            compressed_size += interp_compressed_size;
            return lossless_data;
        }

        
        // compress given the error bound
        uchar *compress_block( Config &conf, T *data, size_t &compressed_size,int tuning=0,int start_level=0,int end_level=0) {
            
            //tuning 0: normal compress 1:tuning to return qbins and psnr 2: tuning to return prediction loss
            /*
            Timer timer;
            timer.start();
            if (conf.trimToZero>0){
                quantizer.setTrimToZero(conf.trimToZero);
            }

            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
            blocksize = conf.interpBlockSize;    
            maxStep=conf.maxStep;    
            interpolator_id = conf.interpMeta.interpAlgo;
            direction_sequence_id = conf.interpMeta.interpDirection;
            alpha=conf.alpha;
            beta=conf.beta;
            size_t cross_block=conf.crossBlock;
            std::vector<uint8_t>interp_ops;
            std::vector<uint8_t>interp_dirs;
            init();
            if (tuning){
                std::vector<int>().swap(quant_inds);
                std::vector<int>().swap(conf.quant_bins);
                conf.quant_bin_counts=std::vector<size_t>(interpolation_level,0);
                conf.decomp_square_error=0.0;
            }
            if(tuning==0 and conf.peTracking){
                prediction_errors.resize(num_elements,0);
                peTracking=1;
            }
            quant_inds.reserve(num_elements);
            size_t interp_compressed_size = 0;
            double eb = quantizer.get_eb();
            if (start_level<=0 or start_level>interpolation_level ){  
                start_level=interpolation_level;     
            } 
            if(end_level>=start_level or end_level<0){
                end_level=0;
            }


            if(!anchor){
                quant_inds.push_back(quantizer.quantize_and_overwrite(*data, 0));
            }
            else if (start_level==interpolation_level){
                if(tuning){
                    conf.quant_bin_counts[start_level-1]=quant_inds.size();
                }
                    
                build_grid(conf,data,maxStep,tuning);
                start_level--;
            }   
            double predict_error=0.0;
            int levelwise_predictor_levels=conf.interpAlgo_list.size();
            for (uint level = start_level; level > end_level && level <= start_level; level--) {
                if (alpha<0) {
                    if (level >= 3) {
                        quantizer.set_eb(eb * eb_ratio);
                    } else {
                        quantizer.set_eb(eb);
                    }
                }
                else if (alpha>=1){           
                    double cur_ratio=pow(alpha,level-1);
                    if (cur_ratio>beta){
                        cur_ratio=beta;
                    }    
                    quantizer.set_eb(eb/cur_ratio);
                }
                else{
                    double cur_ratio=1-(level-1)*alpha;
                    if (cur_ratio<beta){
                        cur_ratio=beta;
                    }            
                    quantizer.set_eb(eb*cur_ratio);
                }
                int cur_interpolator;
                int cur_direction;
                size_t stride = 1U << (level - 1);
                size_t cur_blocksize=blocksize;
                auto inter_block_range = std::make_shared<
                        QoZ::multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                           std::end(global_dimensions),
                                                           cur_blocksize, 0,conf.blockOrder);
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();

                size_t blockwiseSampleBlockSize=(level<=2)?conf.blockwiseSampleBlockSize:cur_blocksize;
                if(blockwiseSampleBlockSize>cur_blocksize)
                    blockwiseSampleBlockSize=cur_blocksize;
                for (auto block = inter_begin; block != inter_end; ++block) {
                    auto start_idx=block.get_global_index();
                    auto end_idx = start_idx;
                    auto sample_end_idx=start_idx;
                    for (int i = 0; i < N; i++) {
                        end_idx[i] += cur_blocksize ;
                        if (end_idx[i] > global_dimensions[i] - 1) {
                            end_idx[i] = global_dimensions[i] - 1;
                        }
                        sample_end_idx[i]+=blockwiseSampleBlockSize;
                        if (sample_end_idx[i] > global_dimensions[i] - 1) {
                            sample_end_idx[i] = global_dimensions[i] - 1;
                        } 
                    }
                    size_t sampled_element_num=1;
                    for(int i=0;i<N;i++){
                        sampled_element_num*=(sample_end_idx[i]-start_idx[i]+1);
                    }
                    std::vector<T> orig_sampled_block(sampled_element_num,0);
                    size_t local_idx=0;
                    if(N==2){
                        for(size_t x=start_idx[0];x<=sample_end_idx[0] ;x++){
                            for(size_t y=start_idx[1];y<=sample_end_idx[1];y++){
                                size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1];
                                orig_sampled_block[local_idx]=data[global_idx];
                                local_idx++;
                            }
                        }
                    }
                    else if(N==3){
                        for(size_t x=start_idx[0];x<=sample_end_idx[0] ;x++){
                            for(size_t y=start_idx[1];y<=sample_end_idx[1] ;y++){
                                for(size_t z=start_idx[2];z<=sample_end_idx[2];z++){
                                    size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1]+z*dimension_offsets[2];
                                    orig_sampled_block[local_idx]=data[global_idx];
                                    local_idx++;
                                }
                            }
                        }
                    }    
                    uint8_t best_op=QoZ::INTERP_ALGO_CUBIC;
                    uint8_t best_dir=0;
                    double best_loss=std::numeric_limits<double>::max();
                    std::vector<int> op_candidates={QoZ::INTERP_ALGO_LINEAR,QoZ::INTERP_ALGO_CUBIC};
                    std::vector<int> dir_candidates={0,QoZ::factorial(N) - 1};
                    for (auto &interp_op:op_candidates) {
                        for (auto &interp_direction: dir_candidates) {
                            double cur_loss=block_interpolation(data, start_idx, sample_end_idx, PB_predict_overwrite,
                                interpolators[interp_op], interp_direction, stride,2);

                            if(cur_loss<best_loss){
                                best_loss=cur_loss;
                                best_op=interp_op;
                                best_dir=interp_direction;
                            }

                            size_t local_idx=0;
                            if(N==2){
                                for(size_t x=start_idx[0];x<=sample_end_idx[0];x++){
                                    for(size_t y=start_idx[1];y<=sample_end_idx[1];y++){
                                        size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1];
                                        data[global_idx]=orig_sampled_block[local_idx];
                                        local_idx++;
                                    }
                                }
                            }
                            else if(N==3){
                                for(size_t x=start_idx[0];x<=sample_end_idx[0];x++){
                                    for(size_t y=start_idx[1];y<=sample_end_idx[1];y++){
                                        for(size_t z=start_idx[2];z<=sample_end_idx[2];z++){
                                            size_t global_idx=x*dimension_offsets[0]+y*dimension_offsets[1]+z*dimension_offsets[2];
                                            data[global_idx]=orig_sampled_block[local_idx];
                                            local_idx++;
                                        }
                                    }
                                }
                            }

                        }
                    }
                    interp_ops.push_back(best_op);
                    interp_dirs.push_back(best_dir);
                    if (peTracking)
                        block_interpolation_block3d(data, start_idx, end_idx, PB_predict_overwrite,
                                    interpolators[best_op], best_dir, stride,3,cross_block);
                    else
                        block_interpolation_block3d(data, start_idx, end_idx, PB_predict_overwrite,
                                    interpolators[best_op], best_dir, stride,0,cross_block);
                }
                if(tuning){
                        
                    conf.quant_bin_counts[level-1]=quant_inds.size();
                }
            }
            //timer.start();
            quantizer.set_eb(eb);
            if(peTracking){
                conf.predictionErrors=prediction_errors;
                conf.interp_ops=interp_ops;
                conf.interp_dirs=interp_dirs;
            }
            if (tuning){
                conf.quant_bins=quant_inds;
                std::vector<int>().swap(quant_inds);
                conf.decomp_square_error=predict_error;
                size_t bufferSize = 1;
                uchar *buffer = new uchar[bufferSize];
                buffer[0]=0;
                
                return buffer;
            }
            if(conf.verbose)
                timer.stop("prediction");  
            //assert(quant_inds.size() == num_elements);
            size_t bufferSize = 2.5 * (quant_inds.size() * sizeof(T) + quantizer.size_est());//ori is 3. In fact this often causes bug.
            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;
            write(global_dimensions.data(), N, buffer_pos);
            write(blocksize, buffer_pos);
            write(interpolator_id, buffer_pos);
            write(direction_sequence_id, buffer_pos);
            write(alpha,buffer_pos);
            write(beta,buffer_pos);         
            write(maxStep,buffer_pos);
            write(levelwise_predictor_levels,buffer_pos);
            write(conf.blockwiseTuning,buffer_pos);
            write(conf.fixBlockSize,buffer_pos);
            write(cross_block,buffer_pos);
            write(conf.trimToZero,buffer_pos);
            write(conf.blockOrder,buffer_pos);
            if(conf.blockwiseTuning){
                size_t ops_num=interp_ops.size();
                write(ops_num,buffer_pos);
                write(interp_ops.data(),ops_num,buffer_pos);
                write(interp_dirs.data(),ops_num,buffer_pos);
            }
            else if(levelwise_predictor_levels>0){
                //write(conf.interpAlgo_list.data(),levelwise_predictor_levels,buffer_pos);
                //write(conf.interpDirection_list.data(),levelwise_predictor_levels,buffer_pos);
                 write(conf.interpMeta_list.data(),levelwise_predictor_levels,buffer_pos);
            }
           
            quantizer.save(buffer_pos);
            quantizer.postcompress_data();
            quantizer.clear();
            encoder.preprocess_encode(quant_inds, 0);
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();
            //timer.stop("Coding");
            //timer.start();
            assert(buffer_pos - buffer < bufferSize);
            uchar *lossless_data = lossless.compress(buffer,
                                                     buffer_pos - buffer,
                                                     compressed_size);
            lossless.postcompress_data(buffer);
            //timer.stop("Lossless") ;
            compressed_size += interp_compressed_size;

            return lossless_data;
            */
            return NULL;
        }
        

        uchar *encoding_lossless(size_t &compressed_size,const std::vector<int> &q_inds=std::vector<int>()) {

            if(q_inds.size()>0)
                quant_inds=q_inds;
            size_t bufferSize = 2.5 * (quant_inds.size() * sizeof(T) + quantizer.size_est());//original is 3
            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;
            quantizer.save(buffer_pos);
            quantizer.clear();
            quantizer.postcompress_data();
            //timer.start();
            encoder.preprocess_encode(quant_inds, 0);
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();
//            timer.stop("Coding");
            assert(buffer_pos - buffer < bufferSize);
            //timer.start();
            uchar *lossless_data = lossless.compress(buffer,
                                                     buffer_pos - buffer,
                                                     compressed_size);
            lossless.postcompress_data(buffer);
//            timer.stop("Lossless");
            return lossless_data;

        }
        void set_eb(double eb){
            quantizer.set_eb(eb);
        }

    private:

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover
        };
        
        void init() {
            assert(blocksize % 2 == 0 && "Interpolation block size should be even numbers");
            num_elements = 1;

            interpolation_level = -1;

            for (int i = 0; i < N; i++) {
                if (interpolation_level < ceil(log2(global_dimensions[i]))) {
                    interpolation_level = (uint) ceil(log2(global_dimensions[i]));
                }
                num_elements *= global_dimensions[i];
            }
            if (maxStep>0){
                anchor=true;//recently moved out of if
                int max_interpolation_level=(uint)log2(maxStep)+1;
                if (max_interpolation_level<=interpolation_level){ 
                    interpolation_level=max_interpolation_level;
                }
            }
            dimension_offsets[N - 1] = 1;
            for (int i = N - 2; i >= 0; i--) {
                dimension_offsets[i] = dimension_offsets[i + 1] * global_dimensions[i + 1];
            }
            dimension_sequences = std::vector<std::array<int, N>>();
            auto sequence = std::array<int, N>();
            for (int i = 0; i < N; i++) {
                sequence[i] = i;
            }
            do {
                dimension_sequences.push_back(sequence);
            } while (std::next_permutation(sequence.begin(), sequence.end()));  
            /*
            mark.clear();
            mark.resize(num_elements,false);
            */
            
        }
       
        void build_grid(Config &conf, T *data,size_t maxStep,int tuning=0){
            
            assert(maxStep>0);

           
            if(tuning>1)
                return;
            /*
            else if(tuning==1 and conf.sampleBlockSize<conf.maxStep and conf.tuningTarget==QoZ::TUNING_TARGET_RD){
                //std::cout<<"dd"<<std::endl;
                quantizer.insert_unpred(*data);
                return;

            }
            */
            if (N==2){
                for (size_t x=maxStep*(tuning==1);x<conf.dims[0];x+=maxStep){
                    for (size_t y=maxStep*(tuning==1);y<conf.dims[1];y+=maxStep){

                        quantizer.insert_unpred(*(data+x*conf.dims[1]+y));
                        /*
                        if(peTracking){
                           // prediction_errors[x*dimension_offsets[0]+y]=*(data+x*dimension_offsets[0]+y);
                            prediction_errors[x*dimension_offsets[0]+y]=0;
                        }*/
                        quant_inds.push_back(0);
                       //mark[x*conf.dims[1]+y]=true;
                    }
                }
            }
            else if(N==3){

                std::array<size_t,3>anchor_strides={maxStep,maxStep,maxStep};

                int fd=conf.fused_dim;
                if(fd>=0)
                    anchor_strides[fd]=1;
                

                
                for (size_t x=anchor_strides[0]*(tuning==1);x<conf.dims[0];x+=anchor_strides[0]){
                    for (size_t y=anchor_strides[1]*(tuning==1);y<conf.dims[1];y+=anchor_strides[1]){
                        for(size_t z=anchor_strides[2]*(tuning==1);z<conf.dims[2];z+=anchor_strides[2]){
                            quantizer.insert_unpred(*(data+x*dimension_offsets[0]+y*dimension_offsets[1]+z) );
                            /*
                            if(peTracking){
                               // prediction_errors[x*dimension_offsets[0]+y*dimension_offsets[1]+z]=*(data+x*dimension_offsets[0]+y*dimension_offsets[1]+z);
                                prediction_errors[x*dimension_offsets[0]+y*dimension_offsets[1]+z]=0;
                            }*/
                            //if(tuning==0)
                                //mark[x*conf.dims[1]*conf.dims[2]+y*conf.dims[2]+z]=true;
                            quant_inds.push_back(0);
                        }           
                    }
                }
            }
        }
 
        void recover_grid(T *decData,const std::array<size_t,N>& global_dimensions,size_t maxStep,size_t fused_dim=-1){
            assert(maxStep>0);
            if (N==2){
                for (size_t x=0;x<global_dimensions[0];x+=maxStep){
                    for (size_t y=0;y<global_dimensions[1];y+=maxStep){
                        decData[x*dimension_offsets[0]+y]=quantizer.recover_unpred();
                        quant_index++;
                    }
                }
            }
            else if(N==3){
                std::array<size_t,3>anchor_strides={maxStep,maxStep,maxStep};
                if(fused_dim>=0)
                    anchor_strides[fused_dim]=1;
                for (size_t x=0;x<global_dimensions[0];x+=anchor_strides[0]){
                    for (size_t y=0;y<global_dimensions[1];y+=anchor_strides[1]){
                        for(size_t z=0;z<global_dimensions[2];z+=anchor_strides[2]){
                            decData[x*dimension_offsets[0]+y*dimension_offsets[1]+z]=quantizer.recover_unpred();
                            quant_index++;
                        }    
                    }
                }

            }
        }
        inline void quantize(size_t idx, T &d, T pred) {
//            preds[idx] = pred;
//            quant_inds[idx] = quantizer.quantize_and_overwrite(d, pred);
            //T orig=d;
            /*
            if(anchor and anchor_threshold>0 and cur_level>=min_anchor_level and fabs(d-pred)>=anchor_threshold){
                quantizer.insert_unpred(d);
                quant_inds.push_back(0);

            }
            else
            */
                quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
            //return fabs(d-orig);
        }

        inline double quantize_tuning(size_t idx, T &d, T pred, int mode=1) {

//            preds[idx] = pred;
//            quant_inds[idx] = quantizer.quantize_and_overwrite(d, pred);

            if (mode==1){
                T orig=d;
                quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
                return (d-orig)*(d-orig);

            }
            else{//} if (mode==2){
                double pred_error=fabs(d-pred);
                /*
                if(peTracking)
                    prediction_errors[idx]=pred_error;
                */
                int q_bin=quantizer.quantize_and_overwrite(d, pred,false);
                /*
                if(peTracking){
                    prediction_errors[idx]=pred_error;
                
                    quant_inds.push_back(q_bin);
                }
                */
                return pred_error;
            }
            /*
            else{
                double pred_error=pred-d;
                int q_bin=quantizer.quantize_and_overwrite(d, pred);
                prediction_errors[idx]=pred_error;
                quant_inds.push_back(q_bin);
                return pred_error;
            }
            */
        }

        inline void recover(size_t idx, T &d, T pred) {
            d = quantizer.recover(pred, quant_inds[quant_index++]);
        };

        inline double quantize_integrated(size_t idx, T &d, T pred, int mode=0){
            double pred_error=0;
            if(mode==-1){//recover
                d = quantizer.recover(pred, quant_inds[quant_index++]);
                return 0;
            }
            else if(mode==0){
                 quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
                 return 0;
            }
            else if(mode==1){
                T orig=d;
                quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
                return (d-orig)*(d-orig);
            }
            else{// if (mode==2){
                pred_error=fabs(d-pred);
                int q_bin=quantizer.quantize_and_overwrite(d, pred,false);
                return pred_error;
            }
        }


        double block_interpolation_1d_cross(T *data, size_t begin, size_t end, size_t stride, const std::string &interp_func, const PredictorBehavior pb,const QoZ::Interp_Meta &meta,int tuning=0,size_t cross_block=0,size_t axis_begin=0,size_t axis_stride=0,size_t cur_axis=0) {//cross block: 0: no cross 1: only front-cross 2: all cross
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }
            double predict_error = 0;
            //uint8_t cubicSplineType=meta.cubicSplineType;

            size_t stride3x = 3 * stride;
            size_t stride5x = 5 * stride;

            if (interp_func == "linear" || (n < 5 and cross_block==0 ) ) {//this place maybe some bug
                if (pb == PB_predict_overwrite) {

                    if (tuning){
                        for (size_t i = 1; i + 1 < n; i += 2) {
                            T *d = data + begin + i * stride;
                            predict_error+=quantize_tuning(d - data, *d, interp_linear(*(d - stride), *(d + stride)),tuning);
                        }
                        if (n % 2 == 0) {
                            size_t offset = begin + (n - 1) * stride;
                            T *d = data + offset;
                            if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis]){
                                predict_error+=quantize_tuning(d - data, *d, interp_linear(*(d - stride), *(d + stride)),tuning);
                            }
                            else if (n >= 4 or (cross_block and  axis_begin >= (4-n)*axis_stride )) {
                                predict_error+=quantize_tuning(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)),tuning);
                            } 
                            else { 
                                predict_error+=quantize_tuning(d - data, *d, *(d - stride),tuning);
                            }
                        }

                    }
                    else{
                        for (size_t i = 1; i + 1 < n; i += 2) {
                            T *d = data + begin + i * stride;
                            quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                        }
                      

                        if (n % 2 == 0) {
                            size_t offset = begin + (n - 1) * stride;
                            T *d = data + offset;
                            if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis]){
                                quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                            }
                            else if (n >= 4 or (cross_block and axis_begin >= (4-n)*axis_stride )) {
                               quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride) ) );                                  
                            } 
                            else { 
                                quantize(d - data, *d, *(d - stride) );
                            }
                        }
                    }
                } else {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        size_t offset = begin + (n - 1) * stride;
                        T *d = data + offset;
                        if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis]){
                            recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)) );
                        }
                        else if (n >= 4 or (cross_block and  axis_begin >= (4-n)*axis_stride )) {
                            recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)) );                                  
                        } 
                        else { 
                            recover(d - data, *d, *(d - stride) );
                        }
                    }
                }
            } else {
                auto interp_cubic=meta.cubicSplineType==0?interp_cubic_1<T>:interp_cubic_2<T>;
                if (pb == PB_predict_overwrite) {

                    if(tuning){
                        T *d;
                        size_t i;
                        for (i = 3; i + 3 < n; i += 2) {
                            d = data + begin + i * stride;
                            predict_error+=quantize_tuning(d - data, *d,
                                     interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),tuning);                           
                        }
                        d = data + begin + stride;
                        if(cross_block and axis_begin >= 2*axis_stride){
                            if(axis_begin+4*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+4*stride<end) )
                                predict_error+=quantize_tuning(d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),tuning);
                            else if (axis_begin+2*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+2*stride<end))
                                predict_error+=quantize_tuning(d - data, *d,
                                    interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride) ),tuning);
                            else
                                predict_error+=quantize_tuning(d - data, *d,
                                    interp_linear1(*(d - stride3x), *(d - stride) ),tuning);                            
                        }
                        else{
                            if(axis_begin+4*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+4*stride<end) )
                                predict_error+=quantize_tuning(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),tuning);
                            else if (axis_begin+2*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+2*stride<end) )
                                predict_error+=quantize_tuning(d - data, *d, interp_linear(*(d - stride), *(d + stride) ),tuning);
                            else
                                predict_error+=quantize_tuning(d - data, *d, *(d - stride) ,tuning);
                        }                       
                        if (begin + i * stride<end){
                            d = data + begin + i * stride;
                            if(cross_block==2 and axis_begin+(i+3)*axis_stride<global_dimensions[cur_axis] and axis_begin+(i-3)*axis_stride>=0)
                                predict_error+=quantize_tuning(d - data, *d,
                                         interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),tuning);
                            else if (axis_begin+(i-3)*axis_stride>=0 and (cross_block>0 or i>=3) ){
                                if (axis_begin+(i+1)*axis_stride<global_dimensions[cur_axis] ){                                   
                                    predict_error+=quantize_tuning(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),tuning );
                                }
                                else if (cross_block>0 or i>=3){
                                    predict_error+=quantize_tuning(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride) ),tuning );
                                }
                                else{
                                    predict_error+=quantize_tuning(d - data, *d, *(d - stride) ,tuning );
                                }
                            }
                            else{
                                if (axis_begin+(i+1)*axis_stride<global_dimensions[cur_axis]){
                                    predict_error+=quantize_tuning(d - data, *d, interp_linear( *(d - stride), *(d + stride)),tuning );
                                }
                                else{
                                    predict_error+=quantize_tuning(d - data, *d, *(d - stride),tuning );
                                }

                            }
                        }             
                        if (n % 2 == 0 and n>4) {
                            size_t offset=begin + (n - 1) * stride;
                            d = data + offset;
                            if(cross_block==2 and axis_begin+(n+2)*axis_stride<global_dimensions[cur_axis] and axis_begin+(n-4)*axis_stride>=0)
                                predict_error+=quantize_tuning(d - data, *d,
                                     interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),tuning);
                            else if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis] and axis_begin+(n-4)*axis_stride>=0)
                                predict_error+=quantize_tuning(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),tuning);
                            else if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis])
                                predict_error+=quantize_tuning(d - data, *d, interp_linear( *(d - stride), *(d + stride)),tuning);
                            else if (axis_begin >= (6-n)*axis_stride and (cross_block>0 or n>=6) )
                                predict_error+=quantize_tuning(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)),tuning);
                            else if (axis_begin >= (4-n)*axis_stride and (cross_block>0 or n>=4) )
                                predict_error+=quantize_tuning(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)),tuning);
                            else
                                predict_error+=quantize_tuning(d - data, *d, *(d - stride),tuning);                          
                        }
                    }                    
                    else{
                        T *d;
                        size_t i;
                        for (i = 3; i + 3 < n; i += 2) {
                            d = data + begin + i * stride;
                            quantize(d - data, *d,
                                     interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)) );
                        }
                        d = data + begin + stride;
                        if(cross_block and axis_begin >= 2*axis_stride){
                            if(axis_begin+4*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+4*stride<end) ){
                                quantize(d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                            }
                            else if (axis_begin+2*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+2*stride<end)) {
                                quantize(d - data, *d,
                                    interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride) ));
                            }
                            else{
                                quantize(d - data, *d,
                                    interp_linear1(*(d - stride3x), *(d - stride) ) );
                            }                          
                        }

                        else{
                            if(axis_begin+4*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+4*stride<end) ){
                                quantize(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)) );
                            }
                            else if (axis_begin+2*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+2*stride<end) ){
                                quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride) ));
                            }
                            else{
                                quantize(d - data, *d, *(d - stride) );
                            }
                        }

                        if (begin + i * stride<end){
                            d = data + begin + i * stride;
                            if(cross_block==2 and axis_begin+(i+3)*axis_stride<global_dimensions[cur_axis] and axis_begin+(i-3)*axis_stride>=0){
                                quantize(d - data, *d,
                                         interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)) );
                            }
                            else if (axis_begin+(i-3)*axis_stride>=0 and (cross_block>0 or i>=3)){
                                if (axis_begin+(i+1)*axis_stride<global_dimensions[cur_axis]){
                                    quantize(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)) );
                                }
                                else if (cross_block>0 or i>=3){                                  
                                    quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride) ) );
                                }
                                else{
                                    quantize(d - data, *d,  *(d - stride)  );
                                }
                            }
                            else{
                                if (axis_begin+(i+1)*axis_stride<global_dimensions[cur_axis]){                                  
                                    quantize(d - data, *d, interp_linear( *(d - stride), *(d + stride)) );
                                }
                                else{
                                    quantize(d - data, *d, *(d - stride) );
                                }
                            }
                        }
                        if (n % 2 == 0 and n>4) {
                            size_t offset=begin + (n - 1) * stride;
                            d = data + offset;
                            if(cross_block==2 and axis_begin+(n+2)*axis_stride<global_dimensions[cur_axis] and axis_begin+(n-4)*axis_stride>=0){
                                quantize(d - data, *d,
                                     interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)) );
                            }
                            else if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis] and axis_begin+(n-4)*axis_stride>=0){
                                quantize(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)) );
                            }
                            else if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis]){
                                quantize(d - data, *d, interp_linear( *(d - stride), *(d + stride)) );
                            }
                            else if (axis_begin >= (6-n)*axis_stride and (cross_block>0 or n>=6)){
                                quantize(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)) );
                            }
                            else if (axis_begin >= (4-n)*axis_stride and (cross_block>0 or n>=4) ){
                                quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)) );
                            }
                            else{
                                quantize(d - data, *d, *(d - stride) );
                            }
                        }                 
                    }
                } else {
                    T *d;
                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        recover(d - data, *d,
                            interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)) );
                    }
                    d = data + begin + stride;
                    if(cross_block and axis_begin >= 2*axis_stride){
                        if(axis_begin+4*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+4*stride<end) ){
                            recover(d - data, *d,
                                interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                        }
                        else if (axis_begin+2*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+2*stride<end)) {
                            recover(d - data, *d,
                                interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride) ));
                        }
                        else{
                             recover(d - data, *d,
                                interp_linear1(*(d - stride3x), *(d - stride) ));
                        }
                    }
                    else{
                        if(axis_begin+4*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+4*stride<end) ){
                            recover(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));
                        }
                        else if(axis_begin+2*axis_stride<global_dimensions[cur_axis] and (cross_block==2 or begin+2*stride<end) ){
                            recover(d - data, *d, interp_linear(*(d - stride), *(d + stride) ));
                        }
                        else{
                            recover(d - data, *d, *(d - stride));
                        }
                    }
                    if (begin + i * stride<end){
                        d = data + begin + i * stride;
                        if(cross_block==2 and axis_begin+(i+3)*axis_stride<global_dimensions[cur_axis] and axis_begin+(i-3)*axis_stride>=0){
                            recover(d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)) );
                        }
                        else if (axis_begin+(i-3)*axis_stride>=0 and (cross_block>0 or i>=3) ){
                            if (axis_begin+(i+1)*axis_stride<global_dimensions[cur_axis]){
                                recover(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)) );
                            }
                            else if (cross_block>0 or i>=3){
                                recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride) ) );
                            }
                            else{
                                recover(d - data, *d,  *(d - stride)  );
                            }
                        }
                        else{
                            if (axis_begin+(i+1)*axis_stride<global_dimensions[cur_axis]){
                                recover(d - data, *d, interp_linear( *(d - stride), *(d + stride)) );
                            }
                            else{
                                recover(d - data, *d, *(d - stride) );
                            }
                        }
                    }
                    if (n % 2 == 0 and n>4) {
                        size_t offset=begin + (n - 1) * stride;
                        d = data + offset;
                        if(cross_block==2 and axis_begin+(n+2)*axis_stride<global_dimensions[cur_axis] and axis_begin+(n-4)*axis_stride>=0){
                            recover(d - data, *d,
                                interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)) );
                        }
                        else if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis] and axis_begin+(n-4)*axis_stride>=0){
                            recover(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)) );
                        }
                        else if (cross_block==2 and axis_begin+n*axis_stride<global_dimensions[cur_axis]){
                            recover(d - data, *d, interp_linear( *(d - stride), *(d + stride)) );

                        }
                        else if (axis_begin >= (6-n)*axis_stride and (cross_block>0 or n>=6) ){
                            recover(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)) );
                        }
                        else if (axis_begin >= (4-n)*axis_stride and (cross_block>0 or n>=4) ){
                            recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)) );
                        }
                        else{
                            recover(d - data, *d, *(d - stride) );
                        }
                    }
                }
            }
            return predict_error;
        }
        
        double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride,const std::string &interp_func,const PredictorBehavior pb,const QoZ::Interp_Meta &meta,int tuning=0) {
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }
            double predict_error = 0;
            size_t stride2x=2*stride;
            size_t stride3x = 3 * stride;
            size_t stride5x = 5 * stride;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            if (interp_func == "linear" || n < 5) {
                for (size_t i = 1; i + 1 < n; i += 2) {
                    T *d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if (n < 4) {                              
                        predict_error+=quantize_integrated(d - data, *d, *(d - stride),mode);
                        } else {
                        predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)),mode);
                    }
                }
            } else {
                auto interp_cubic=meta.cubicSplineType==0?interp_cubic_1<T>:interp_cubic_2<T>;
                T *d;
                size_t i;
                if(!meta.adjInterp){
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        predict_error+=quantize_integrated(d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                    }

                    d = data + begin + stride;

                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),mode);
                    d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        predict_error+=
                        quantize_integrated(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)),mode);
                    }
                }
                else{// if(meta.adjInterp==1){
                    auto interp_cubic_adj=meta.cubicSplineType==0?interp_cubic_adj_2<T>:interp_cubic_adj_1<T>;
                    //i=1
                    d = data + begin + stride;
                    
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),mode);
                    //predict_error+=quantize_integrated(d - data, *d, interp_cubic_front_adj(*(d -stride),*(d + stride), *(d+stride2x), *(d + stride3x)),mode);
                    for (i = 5; i + 3 < n; i += 4) {
                        
                        d = data + begin + i * stride;
                     
                        predict_error+=quantize_integrated(d - data, *d,
                                    interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                    }

                    //i=n-3 or n-2

                    if(i<n-1){
                        d = data + begin + i * stride;
                       
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                        //predict_error+=quantize_integrated(d - data, *d, interp_cubic_back_adj(*(d -stride3x),*(d - stride2x), *(d-stride), *(d + stride)),mode);
                    }
                    //i=n-1
                    else if(i<n){
                         d = data + begin + (n - 1) * stride;
             
                        predict_error+=
                        //quantize_integrated(d - data, *d, interp_quad_3_adj(*(d - stride3x), *(d - stride2x), *(d - stride)),mode);
                        quantize_integrated(d - data, *d, *(d - stride),mode);//to determine
                        //quantize_integrated(d - data, *d, *(d - stride),mode);

                    }


                    for (i = 3; i + 3 < n; i += 4) {
                        d = data + begin + i * stride;

                        predict_error+=quantize_integrated(d - data, *d,
                                   interp_cubic_adj(*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x),*(d + stride3x)),mode);
                        //predict_error+=quantize_integrated(d - data, *d,
                         //           interp_cubic_3(*(d - stride2x), *(d - stride), *(d + stride), *(d+stride2x)),mode);
                    }
                    //i=n-3 or n-2
                    if(i<n-1){
                        d = data + begin + i * stride;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x), *(d - stride), *(d + stride)),mode);
                        //predict_error+=quantize_integrated(d - data, *d, interp_cubic_back_adj(*(d -stride3x),*(d - stride2x), *(d-stride), *(d + stride)),mode);
                    }
                    //i=n-1
                    else if(i<n){
                         d = data + begin + (n - 1) * stride;

                        predict_error+=
                        //quantize_integrated(d - data, *d, interp_quad_3_adj(*(d - stride3x), *(d - stride2x), *(d - stride)),mode);
                        quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x),*(d - stride)) ,mode);//to determine
                        //quantize_integrated(d - data, *d, *(d - stride),mode);
                    }
                }
                /*
                else if(meta.adjInterp==2){
                    auto interp_cubic_adj=meta.cubicSplineType==0?interp_cubic_adj_4<T>:interp_cubic_adj_3<T>;
                    d = data + begin + stride;

                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),mode);
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        predict_error+=quantize_integrated(d - data, *d,
                                    interp_cubic_adj(*(d - stride3x),*(d - stride2x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                    }

                    
                    d = data + begin + i * stride;
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x), *(d - stride), *(d + stride)),mode);
                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        predict_error+=
                        quantize_integrated(d - data, *d, lorenzo_1d(*(d - stride2x), *(d - stride)),mode);
                    }

                }
                */

                
                
                
                
            }
            return predict_error;
        } 


        double regressive_interpolation_1d_linear(T *data, int &status, const size_t & cur_idx,const size_t &main_direction,const std::vector<size_t> & sub_directions, const std::vector<size_t> & strides, const std::vector<size_t> & strides3x,const std::vector<size_t> &dimensional_sparsity){
            T *d=data+cur_idx;
            std::vector<double> A={*(d-strides3x[main_direction]),*(d+strides[main_direction]),*(d-strides[main_direction]),*(d+strides3x[main_direction])};
            std::vector<double> b={*(d-strides[main_direction]),*(d+strides[main_direction])};
            /*
            for(auto sub_direction:sub_directions){
                size_t sub_stride=dimensional_sparsity[sub_direction]*strides[sub_direction];
                std::vector<double>tempA={*(d-sub_stride-strides3x[main_direction]),*(d-sub_stride+strides[main_direction]),*(d-sub_stride-strides[main_direction]),*(d-sub_stride+strides3x[main_direction])
                    ,*(d+sub_stride-strides3x[main_direction]),*(d+sub_stride+strides[main_direction]),*(d+sub_stride-strides[main_direction]),*(d+sub_stride+strides3x[main_direction])};
                A.insert(A.end(),tempA.begin(),tempA.end());
                

                std::vector<double>tempb={*(d-sub_stride-strides[main_direction]),*(d-sub_stride+strides[main_direction]),*(d+sub_stride-strides[main_direction]),*(d+sub_stride+strides[main_direction])};
                b.insert(b.end(),tempb.begin(),tempb.end());

            }
        */
            auto reg_res=QoZ::Regression(A.data(),b.size(),2,b.data(),status);
            if(status==0){
                if(isnan(reg_res[0]) or isnan(reg_res[1]) or fabs(reg_res[0])>0.75 or fabs(reg_res[1])>0.75 or fabs(reg_res[0])<0 or fabs(reg_res[1])<0){
                    status=1;
                    return 0;
                }
                else
                    return reg_res[0]*(*(d-strides[main_direction]))+reg_res[1]*(*(d+strides[main_direction]));
            }
            else{
                return 0;
            }
        }

        double regressive_interpolation_1d_cubic(T *data, int &status, const size_t & cur_idx,const size_t &main_direction,const std::vector<size_t> & sub_directions, const std::vector<size_t> & strides, const std::vector<size_t> & strides3x,const std::vector<size_t> &dimensional_sparsity){
            T *d=data+cur_idx;
            size_t main_stride=strides[main_direction],main_stride3x=main_stride*3,main_stride5x=main_stride*5,main_stride7x=main_stride*7,main_stride9x=main_stride*9;
            std::vector<double> A={*(d-main_stride9x),*(d-main_stride5x),*(d-main_stride),*(d+main_stride3x),
                *(d-main_stride7x),*(d-main_stride3x),*(d+main_stride),*(d+main_stride5x),
                *(d-main_stride5x),*(d-main_stride),*(d+main_stride3x),*(d+main_stride7x),
                *(d-main_stride3x),*(d+main_stride),*(d+main_stride3x),*(d+main_stride7x)};
            std::vector<double> b={*(d-main_stride3x),*(d-main_stride),*(d+main_stride),*(d+main_stride3x)};
            /*
            for(auto sub_direction:sub_directions){
                size_t sub_stride=dimensional_sparsity[sub_direction]*strides[sub_direction];
                std::vector<double>tempA={*(d-sub_stride-strides3x[main_direction]),*(d-sub_stride+strides[main_direction]),*(d-sub_stride-strides[main_direction]),*(d-sub_stride+strides3x[main_direction])
                    ,*(d+sub_stride-strides3x[main_direction]),*(d+sub_stride+strides[main_direction]),*(d+sub_stride-strides[main_direction]),*(d+sub_stride+strides3x[main_direction])};
                A.insert(A.end(),tempA.begin(),tempA.end());
                

                std::vector<double>tempb={*(d-sub_stride-strides[main_direction]),*(d-sub_stride+strides[main_direction]),*(d+sub_stride-strides[main_direction]),*(d+sub_stride+strides[main_direction])};
                b.insert(b.end(),tempb.begin(),tempb.end());

            }
            */
            auto reg_res=QoZ::Regression(A.data(),b.size(),4,b.data(),status);
            if(status==0){
                
                for(size_t i=0;i<4;i++){
                    if(isnan(reg_res[i]) or fabs(reg_res[i])>0.9 or fabs(reg_res[i])<-0.4){
                        status=1;
                        return 0;
                    }
                }
                return reg_res[0]*(*(d-main_stride3x))+reg_res[1]*(*(d-main_stride))+reg_res[2]*(*(d+main_stride))+reg_res[3]*(*(d+main_stride3x));
            }
            else{
                return 0;
            }
        }
        double block_interpolation_1d_regressive(T *data, const std::vector<size_t> &block_begins,const std::vector<size_t> & block_ends, const size_t & main_direction,const std::vector<size_t> &begins,const std::vector<size_t> & ends,const std::vector<size_t> &dimensional_sparsity,size_t m_stride,const std::string &interp_func,const PredictorBehavior pb,const QoZ::Interp_Meta &meta,int tuning=0) {
            //all coords and strides are mathematical.
            //ends are included in the interpolations.
            //dimensional offsets are in the global scale.
            //dimensional_sparsity: {2,2,2} to {1,1,1}
            size_t num_dims=block_begins.size();
            std::vector<size_t> strides(num_dims,0);
            std::vector<size_t> strides3x(num_dims,0);
            std::vector<size_t> strides5x(num_dims,0);
            for(size_t i=0;i<N;i++){
                strides[i]=m_stride*dimension_offsets[i];
                strides3x[i]=strides[i]*3;
                strides5x[i]=strides[i]*3;
            }
            size_t stride=strides[main_direction],stride3x=strides3x[main_direction],stride5x=strides5x[main_direction];
            std::vector<size_t> sub_directions;
            for(size_t i=0;i<N;i++){
                if(i!=main_direction)
                    sub_directions.push_back(i);
            }

            size_t n = (ends[main_direction] - begins[main_direction])/m_stride + 1;
            if (n <= 1) {
                return 0;
            }

            double predict_error = 0;
            size_t begin=0;
            for(size_t i=0;i<N;i++){
                begin+=begins[i]*dimension_offsets[i];
            }

            bool reg_along_sub_d=true;
            //std::cout<<block_begins[0]<<" "<<block_begins[1]<<" "<<block_ends[0]<<" "<<block_ends[1]<<std::endl;
            for(auto i:sub_directions){
                if (begins[i]<block_begins[i]+dimensional_sparsity[i]*m_stride or ends[i]>block_ends[i]-dimensional_sparsity[i]*m_stride){
                    reg_along_sub_d=false;
                    break;
                }
            }

            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            if (interp_func == "linear" || n < 5) {
               
                for (size_t i = 1; i + 1 < n; i += 2) {

                    size_t cur_idx=begin + i * stride;
                    T *d = data + cur_idx;

                    if(reg_along_sub_d and i>=3 and i+3<n){
                        //std::cout<<begins[0]<<" "<<begins[1]<<" "<<i<<std::endl;
                        //std::vector<size_t> cur_coord=begins;
                        //cur_coord[main_direction]+=i;
                        int status=0;
                        T prediction=regressive_interpolation_1d_linear(data,status,cur_idx,main_direction,sub_directions,strides,strides3x,dimensional_sparsity);
                        /*
                        if(status!=0){
                            std::cout<<begins[0]<<" "<<begins[1]<<" "<<i<<std::endl; 
                        }
                        */
                        //if (prediction<0 or prediction>1)
                            //std::cout<<cur_idx<<" "<<prediction<<" "<<status<<std::endl;
                            
                        //prediction=0;
                            

                        if(status!=0 )
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                        else
                            predict_error+=quantize_integrated(d - data, *d, prediction,mode);

                    }

                    else
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride), *(d + stride)),mode);
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if (n < 4) {                              
                        predict_error+=quantize_integrated(d - data, *d, *(d - stride),mode);
                        } else {
                        predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)),mode);
                    }
                }
            } else {
                auto interp_cubic=meta.cubicSplineType==0?interp_cubic_1<T>:interp_cubic_2<T>;

                T *d;
                size_t i;
                for (i = 3; i + 3 < n; i += 2) {


                    size_t cur_idx=begin + i * stride;
                    d = data + cur_idx;

                    if(reg_along_sub_d and i>=9 and i+9<n){
                        //std::cout<<begins[0]<<" "<<begins[1]<<" "<<i<<std::endl;
                        //std::vector<size_t> cur_coord=begins;
                        //cur_coord[main_direction]+=i;
                        int status=0;
                        T prediction=regressive_interpolation_1d_cubic(data,status,cur_idx,main_direction,sub_directions,strides,strides3x,dimensional_sparsity);
                        if(status!=0 )
                            predict_error+=quantize_integrated(d - data, *d,
                                interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                        else
                            predict_error+=quantize_integrated(d - data, *d, prediction,mode);
                    }
                    else
                        predict_error+=quantize_integrated(d - data, *d,
                                interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)),mode);
                }
                d = data + begin + stride;
                predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)),mode);
                d = data + begin + i * stride;
                predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)),mode);
                if (n % 2 == 0) {
                    d = data + begin + (n - 1) * stride;
                    predict_error+=
                    quantize_integrated(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)),mode);
                }
                
            }
            return predict_error;
        } 
        double block_interpolation_2d(T *data, size_t begin1, size_t end1, size_t begin2, size_t end2, size_t stride1,size_t stride2,const std::string &interp_func,const PredictorBehavior pb,const std::array<double,2> &dim_coeffs,const QoZ::Interp_Meta &meta,int tuning=0) {
            size_t n = (end1 - begin1) / stride1 + 1;
            if (n <= 1) {
                return 0;
            }
            size_t m = (end2 - begin2) / stride2 + 1;
            if (m <= 1) {
                return 0;
            }
            double predict_error = 0;
            
            double coeff_x=(dim_coeffs[0])/((dim_coeffs[0])+(dim_coeffs[1])),coeff_y=1-coeff_x;
            //std::cout<<coeff_x<<" "<<coeff_y<<std::endl;
            //coeff_x=0.5; coeff_y=0.5;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            if (interp_func == "linear"||n<5 ||m<5) {//nm cond temp added
               
        
                for (size_t i = 1; i + 1 < n; i += 2) {
                    for(size_t j=1;j+1<m;j+=2){
                        T *d = data + begin1 + i* stride1+begin2+j*stride2;
                        //std::cout<<"q1 "<<i<<" "<<j<<std::endl;
                        //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1), *(d + stride1),*(d - stride2), *(d + stride2)),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_linear(*(d - stride1), *(d + stride1))+coeff_y*interp_linear(*(d - stride2), *(d + stride2)),mode);

                    }
                    if(m%2 ==0){
                        // std::cout<<"q2 "<<i<<std::endl;
                        T *d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1), *(d + stride1)),mode);//to determine whether 2d or 1d 
                    }
                }
                if (n % 2 == 0) {
                    for(size_t j=1;j+1<m;j+=2){
                        // std::cout<<"q3 "<<j<<std::endl;
                        T *d = data + begin1 + (n-1) * stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride2), *(d + stride2)),mode);//to determine whether 2d or 1d 
                    }
                    if(m%2 ==0){
                        //std::cout<<"q4"<<std::endl;
                        T *d = data + begin1 + (n-1) * stride1+begin2+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d - stride1-stride2), *(d - stride1), *(d - stride2)),mode);//to determine whether use lorenzo or not
                    }          
                }
                    
            }
                    
            else{//cubic
                auto interp_cubic=meta.cubicSplineType==0?interp_cubic_1<T>:interp_cubic_2<T>;
                size_t stride3x1=3*stride1,stride3x2=3*stride2,stride5x1=5*stride1,stride5x2=5*stride2,stride2x1=2*stride1,stride2x2=2*stride2;
                //adaptive todo
              
                   
                size_t i,j;
                T *d;
                if(!meta.adjInterp){
                    for (i = 3; i + 3 < n; i += 2) {
                       
                        for(j=3;j+3<m;j+=2){
                            d = data + begin1 + i* stride1+begin2+j*stride2;


                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                        ,interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //j=1
                        d = data + begin1 + i* stride1+ begin2+stride2;
                        //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                        //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);//to determine
                        //j=m-3 or m-2
                        d = data +begin1 + i* stride1+ begin2+j*stride2;
                        //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                        //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d,interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        //j=m-1
                        if(m%2 ==0){
                            d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }
                    }
                    //i=1
                    for(j=3;j+3<m;j+=2){
                        d = data + begin1 + stride1+begin2+j*stride2;
                        
                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d,interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=1
                    d = data + begin1 + stride1+ begin2+stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                            
                    //j=m-3 or m-2
                    d = data +begin1 + stride1+ begin2+j*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                    //j=m-1
                    if(m%2 ==0){
                        d = data + begin1 + stride1+begin2+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        d = data + begin1 + i*stride1+begin2+j*stride2;
                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)) ,interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);


                    }
                    //j=1
                    d = data + begin1 + i*stride1+ begin2+stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                    
                    //j=m-3 or m-2
                    d = data +begin1 + i*stride1+ begin2+j*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                    //j=m-1
                    if(m%2 ==0){
                        d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)) ,interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=1
                        d = data + begin1 + (n-1)*stride1+ begin2+stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        //j=m-3 or m-2
                        d = data +begin1 + (n-1)*stride1+ begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        //j=m-1
                        if(m%2 ==0){
                            d = data + begin1 + (n-1) * stride1+begin2+(m-1)*stride2;
                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)), interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)) ,mode);
                        } 
                    }
                }
                else{
                    auto interp_cubic_adj=meta.cubicSplineType==0?interp_cubic_adj_2<T>:interp_cubic_adj_1<T>;
                    size_t j_start;
                    //first half (non-adj)
                    //std::cout<<"f1"<<std::endl;
                    for (i = 3; i + 3 < n; i += 2) {
                        j_start= (i%4==1)?5:3;
                        for(j=j_start;j+3<m;j+=4){

 
                            d = data + begin1 + i* stride1+begin2+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                        ,interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d , coeff_x*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=1
                        if(j_start==5){
                            
                            d = data + begin1 + i* stride1+ begin2+stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d,interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);//to determine
                        }
                        
                        //j=m-3 or m-2 or j=m-1
                        if(j<m){
                            d = data +begin1 + i* stride1+ begin2+j*stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                            //predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                ,mode);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);

                        }
                        
                            /*                          
                        //j=m-1
                        else if(j<m){
                            d = data + begin1 + i * stride1+begin2+j*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }

                        */
                        
                    }
                    //std::cout<<"f2"<<std::endl;
                    //i=1
                    for(j=5;j+3<m;j+=4){
                        d = data + begin1 + stride1+begin2+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=1
                    d = data + begin1 + stride1+ begin2+stride2;

                    //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                    //j=m-3 or m-2
                    if(j<m-1){
                        d = data +begin1 + stride1+ begin2+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                    }
                    else if(j<m){//j=m-1

                        d = data + begin1 + stride1+begin2+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    }


                    //i=n-3 or n-2
                     //std::cout<<"f3"<<std::endl;
                    j_start= (i%4==1)?5:3;
                    for(j=j_start;j+3<m;j+=4){
   
                        d = data + begin1 + i*stride1+begin2+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)) ,interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);

                    }
                    
                    //j=1
                    if(j_start==5){
    
                        d = data + begin1 + i*stride1+ begin2+stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                    }
                    //j=m-3 or m-2
                    if(j<m-1){
  
                        d = data +begin1 + i*stride1+ begin2+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))+coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                    }                    
                    //j=m-1
                    else if(j<m){
   
                        d = data + begin1 + i * stride1+begin2+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                    }


                    //i=n-1 (odd)
                    // std::cout<<"f4"<<std::endl;
                    if (n % 2 == 0) {
                        j_start= ((n-1)%4==1)?5:3;
                        for(j=j_start;j+3<m;j+=4){
 
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)) ,interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }

                        //j=1
                        if(j_start==5){
 
                            d = data + begin1 + (n-1)*stride1+ begin2+stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
 
                            d = data +begin1 + (n-1)*stride1+ begin2+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        else if(j<m){
 
                            d = data + begin1 + (n-1) * stride1+begin2+j*stride2;

                            //redict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)), interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),mode);
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)), mode);
                        } 
                    }

                    //second half (adj)
                    // std::cout<<"f5"<<std::endl;
                    for (i = 3; i + 3 < n; i += 2) {
                        j_start= (i%4==1)?3:5;
                        for(j=j_start;j+3<m;j+=4){
                            
                            d = data + begin1 + i* stride1+begin2+j*stride2;

                            //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1))
                            //                                        ,interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  );,mode);
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1))
                                                                            +coeff_y*interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2)) ,mode);
                        }
                        //j=1
                        /*
                        if(mode==-1){
                                std::cout<<i<<" "<<1<<std::endl;
                            }
                        */
                        if(j_start==5){
                     
                            d = data + begin1 + i* stride1+ begin2+stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1)),mode);//to determine
                        }
                        /*
                        if(mode==-1){
                                std::cout<<i<<" "<<j<<std::endl;
                            }
                            */
                        //j=m-3 or m-2 or m-1
                        if(j<m){
         
                            d = data +begin1 + i* stride1+ begin2+j*stride2;

                            //predict_error+=quantize_tuning(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                                , interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1),*(d + stride2x1), *(d + stride3x1)),mode);//to determine
                        }
                        /*
                        //j=m-1
                        else if(j<m){
                            d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                            //predict_error+=quantize_tuning(d - data, *d, interp_linear(interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                            //                     , interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)) ),tuning);
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                         ,mode);
                        }
                        */

                    }

                    //i=1
                    // std::cout<<"f6"<<std::endl;
                    for(j=3;j+3<m;j+=4){
                      
                        d = data + begin1 + stride1+begin2+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2)) ,mode);
                    }
                    /*
                    //j=1
                    d = data + begin1 + stride1+ begin2+stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ),mode);//bug when no nmcond and n or m<=4, all the following quads has this problem
                    */    
                    //j=m-3 or m-2
                    if(j<m-1){
                       
                        d = data +begin1 + stride1+ begin2+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)), interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    //j=m-1
                    else if(j<m){
                        
                        d = data + begin1 + stride1+begin2+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),mode);//to determine
                    }

                    //i= n-3 or n-2
                    // std::cout<<"f7"<<std::endl;
                    j_start= (i%4==1)?3:5;
                    for(j=j_start;j+3<m;j+=4){
                        
                        d = data + begin1 + i* stride1+begin2+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  ,mode);
                    }
                    //j=1
                    if(j_start==5){
                       
                        d = data + begin1 + i*stride1+ begin2+stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                        //                                                            , interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)),mode);
                    }
                    //j=m-3 or m-2
                    if(j<m-1){
                       
                        d = data +begin1 + i*stride1+ begin2+j*stride2;

                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)), interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    //j=m-1
                    else if(j<m){
                       
                        d = data + begin1 + i * stride1+begin2+j*stride2;

                        predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)) ,mode);//to determine
                    }
                    
                    //i==n-1
                    // std::cout<<"f8"<<std::endl;
                    if (n % 2 == 0) {
                        j_start= ((n-1)%4==1)?3:5;
                        for(j=j_start;j+3<m;j+=4){
                            
                            d = data + begin1 + (n-1)* stride1+begin2+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2),*(d + stride2x2), *(d + stride3x2))  ,mode);
                        }
                        //j=1
                        if(j_start==5){
                            
                            d = data + begin1 + (n-1)*stride1+ begin2+stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        if(j<m-1){
                            d = data +begin1 + (n-1)*stride1+ begin2+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ,mode);
                        }
                        //j=m-1
                        else if(j<m){
                            d = data + begin1 + (n-1) * stride1+begin2+j*stride2;

                            predict_error+=quantize_integrated(d - data, *d, lorenzo_2d(*(d-stride1-stride2),*(d-stride1),*(d-stride2)),mode);
                        } 
                    }
                }  
            }      
            return predict_error;
        }
        double block_interpolation_3d(T *data, size_t begin1, size_t end1, size_t begin2, size_t end2, size_t begin3, size_t end3, size_t stride1,size_t stride2,size_t stride3, const std::string &interp_func, const PredictorBehavior pb,const std::array<double,3> &dim_coeffs,const QoZ::Interp_Meta &meta,int tuning=0) {
            size_t n = (end1 - begin1) / stride1 + 1;
            if (n <= 1) {
                return 0;
            }
            size_t m = (end2 - begin2) / stride2 + 1;
            if (m <= 1) {
                return 0;
            }
            size_t p = (end3 - begin3) / stride3 + 1;
            if (p <= 1) {
                return 0;
            }
            double predict_error = 0;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;

            double coeff_x=dim_coeffs[0]/(dim_coeffs[0]+dim_coeffs[1]+dim_coeffs[2]);
            double coeff_y=dim_coeffs[1]/(dim_coeffs[0]+dim_coeffs[1]+dim_coeffs[2]);
            double coeff_z=1-coeff_x-coeff_y;

            double coeff_x_xy=(coeff_x)/(coeff_x+coeff_y),coeff_y_xy=1-coeff_x_xy;
            double coeff_x_xz=(coeff_x)/(coeff_x+coeff_z),coeff_z_xz=1-coeff_x_xz;
            double coeff_y_yz=(coeff_y)/(coeff_y+coeff_z),coeff_z_yz=1-coeff_y_yz;


            if (interp_func == "linear" || n<5 || m<5 || p<5 ){//nmpcond temp added
                
                for (size_t i = 1; i + 1 < n; i += 2) {
                    for(size_t j=1;j+1<m;j+=2){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_linear(*(d - stride1), *(d + stride1))
                                                                            +coeff_y*interp_linear(*(d - stride2), *(d + stride2))
                                                                            +coeff_z*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin1 + i* stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_linear(*(d - stride1), *(d + stride1))
                                                                                +coeff_y_xy*interp_linear(*(d - stride2), *(d + stride2)),mode);
                        }

                    }
                    if(m%2 ==0){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin1 + i* stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_linear(*(d - stride1), *(d + stride1))
                                                                                +coeff_z_xz*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin1 + i* stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1), *(d + stride1)),mode);
                        }
                    }      
                }
                if (n % 2 == 0) {
                    for(size_t j=1;j+1<m;j+=2){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin1 + (n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_linear(*(d - stride2), *(d + stride2))
                                                                                +coeff_z_yz*interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin1 + (n-1)* stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride2), *(d + stride2)),mode);
                        }
                    }
                    if(m%2 ==0){
                        for(size_t k=1;k+1<p;k+=2){
                            T *d = data + begin1 + (n-1)* stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride3), *(d + stride3)),mode);
                        }
                        if(p%2==0){
                            T *d = data + begin1 + (n-1)* stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                        }
                    }           
                }
            }
            else{//cubic
                auto interp_cubic=meta.cubicSplineType==0?interp_cubic_1<T>:interp_cubic_2<T>;
                size_t stride3x1=3*stride1,stride3x2=3*stride2,stride5x1=5*stride1,stride5x2=5*stride2,stride3x3=3*stride3,stride5x3=5*stride3,stride2x1=2*stride1,stride2x2=2*stride2,stride2x3=2*stride3;
                //adaptive todo
              
                   
                size_t i,j,k;
                T *d;
                if(!meta.adjInterp){
                    for (i = 3; i + 3 < n; i += 2) {
                        for(j=3;j+3<m;j+=2){
                            for(k=3;k+3<p;k+=2){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                            }
                            //k=1
                            d = data + begin1 + i* stride1+begin2+j*stride2+begin3+stride3;
                            /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);//should or how we ave for different interps?
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y_xy*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        
                            //k=p-3 or p-2
                            d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y_xy*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            
                            //k=p-1
                            if(p%2==0){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                                */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y_xy*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                                 
                            }
                        }
                        //j=1
                        for(k=3;k+3<p;k+=2){
                            d = data + begin1 + i* stride1+ begin2+stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                         }
                        //k=1
                        d = data + begin1 + i* stride1+begin2+stride2+begin3+stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);

                        //k=p-3 or p-2
                        d = data + begin1 + i* stride1+begin2+stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin1 + i* stride1+begin2+stride2+begin3+(p-1)*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                        //j=m-3 or m-2
                        for(k=3;k+3<p;k+=2){
                            d = data + begin1 + i* stride1+ begin2+j*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                         }
                        //k=1
                        d = data + begin1 + i* stride1+begin2+j*stride2+begin3+stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        //k=p-3 or p-2
                        d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_tuning(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin1 + i* stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                            /*
                            predict_error+=quantize_tuning(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            for(k=3;k+3<p;k+=2){
                                d = data + begin1 + i* stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                    interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                    interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                                */
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1

                            d = data + begin1 + i* stride1+begin2+(m-1)*stride2+begin3+stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                    interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            //k=p-3 or p-2
                            d = data + begin1 + i* stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                    interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            //k=p-1
                            if(p%2==0){
                                d = data + begin1 + i* stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),
                                    interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                    interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                                */
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                        }
                    }
                    //i=1
                    for(j=3;j+3<m;j+=2){
                        for(k=3;k+3<p;k+=2){
                            d = data + begin1 +  stride1+begin2+j*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin1 +  stride1+begin2+j*stride2+begin3+stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        //k=p-3 or p-2
                        d = data + begin1 + stride1+begin2+j*stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin1 +  stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                                interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    for(k=3;k+3<p;k+=2){
                        d = data + begin1 +  stride1+ begin2+stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                            interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                            interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    d = data + begin1 + stride1+begin2+stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                        +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                    //k=p-3 or p-2
                    d = data + begin1 +  stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                        +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    //k=p-1
                    if(p%2==0){
                        d = data + begin1 +  stride1+begin2+stride2+begin3+(p-1)*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=m-3 or m-2
                    for(k=3;k+3<p;k+=2){
                        d = data + begin1 +  stride1+ begin2+j*stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), 
                            interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),
                            interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    d = data + begin1 + stride1+begin2+j*stride2+begin3+stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                    +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    //k=p-3 or p-2
                    d = data + begin1 +  stride1+begin2+j*stride2+begin3+k*stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                    +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                    +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                    //k=p-1
                    if(p%2==0){
                        d = data + begin1 + stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        for(k=3;k+3<p;k+=2){
                            d = data + begin1 +  stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin1 +  stride1+begin2+(m-1)*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        //k=p-3 or p-2
                        d = data + begin1 + stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin1 + stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        for(k=3;k+3<p;k+=2){
                            d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                                interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=1
                        d = data + begin1 +  i*stride1+begin2+j*stride2+begin3+stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                                interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        //k=p-3 or p-2
                        d = data + begin1 +i* stride1+begin2+j*stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                                interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin1 +  i*stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                                interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    for(k=3;k+3<p;k+=2){
                        d = data + begin1 +  i*stride1+ begin2+stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),
                            interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                            interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    d = data + begin1 + i*stride1+begin2+stride2+begin3+stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                    +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                    +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    //k=p-3 or p-2
                    d = data + begin1 +  i*stride1+begin2+stride2+begin3+k*stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                    +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                    +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    //k=p-1
                    if(p%2==0){
                        d = data + begin1 +  i*stride1+begin2+stride2+begin3+(p-1)*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=m-3 or m-2
                    for(k=3;k+3<p;k+=2){
                        d = data + begin1 +  i*stride1+ begin2+j*stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), 
                            interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),
                            interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    d = data + begin1 + i*stride1+begin2+j*stride2+begin3+stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                    +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                    +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                    //k=p-3 or p-2
                    d = data + begin1 +  i*stride1+begin2+j*stride2+begin3+k*stride3;
                    predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                    +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                    +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    //k=p-1
                    if(p%2==0){
                        d = data + begin1 + i*stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        for(k=3;k+3<p;k+=2){
                            d = data + begin1 +  i*stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)), 
                                interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin1 +  i*stride1+begin2+(m-1)*stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        //k=p-3 or p-2
                        d = data + begin1 + i*stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin1 + i*stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            for(k=3;k+3<p;k+=2){
                                d = data + begin1 + (n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                    interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                                */
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                                +coeff_z_yz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                    interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            //k=p-3 or p-2
                            d = data + begin1 +(n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                    interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            //k=p-1
                            if(p%2==0){
                                d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                                /*
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                    interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                                    interp_quad_3(*(d - stride5x3), *(d - stride3x3), *(d - stride3)) ),mode);
                                */
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                        }
                        //j=1
                        for(k=3;k+3<p;k+=2){
                            d = data + begin1 + (n-1)*stride1+ begin2+stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin1 + (n-1)*stride1+begin2+stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                        +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        //k=p-3 or p-2
                        d = data + begin1 +  (n-1)*stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                        +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin1 + (n-1)*stride1+begin2+stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //j=m-3 or m-2
                        for(k=3;k+3<p;k+=2){
                            d = data + begin1 + (n-1)*stride1+ begin2+j*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_3(*(d - stride5x1), *(d - stride3x1), *(d - stride1)),
                                interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        d = data + begin1 + (n-1)*stride1+begin2+j*stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))
                                                                        +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        //k=p-3 or p-2
                        d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                        +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                        //k=p-1
                        if(p%2==0){
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            for(k=3;k+3<p;k+=2){
                                d = data + begin1 +  (n-1)*stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            d = data + begin1 +  (n-1)*stride1+begin2+(m-1)*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);

                            //k=p-3 or p-2
                            d = data + begin1 + (n-1)*stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            //k=p-1
                            if(p%2==0){
                                d = data + begin1 + (n-1)*stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }
                }
                else{
                    auto interp_cubic_adj=meta.cubicSplineType==0?interp_cubic_adj_2<T>:interp_cubic_adj_1<T>;

                    size_t k_start;
                    //first half (non-adj) 

                    //for(size_t round=0;round<=1;round++){

                    size_t ks1=3,ks2=5;
                    /*
                    auto interp_cubic_func=round==0?interp_cubic_noadj:interp_cubic_adj;
                    auto interp_quad_1_func=round==0?interp_quad_1<T>:interp_quad_1_adj<T>;
                    auto interp_quad_2_func=round==0?interp_quad_2<T>:interp_quad_2_adj<T>;
                    */


                    for (i = 3; i + 3 < n; i += 2) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                            }
                        
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                            }

                        }
                        //j=1
                        k_start=(i+1)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + i* stride1+ begin2+stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + i* stride1+begin2+stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }

                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin1 + i* stride1+begin2+stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }

                        //j=m-3 or m-2 or m-1
                        while(j<m){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + i* stride1+ begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin1 + i* stride1+begin2+j*stride2+begin3+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                            }
                            j+=2;
                        }
                    }
                    //i=1
                    
                    for(j=3;j+3<m;j+=2){
                        k_start=(1+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  stride1+begin2+j*stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  stride1+begin2+j*stride2+begin3+stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if (k<p){
                            d = data + begin1 + stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
    
                    }
                    //j=1 (i+j=2)
                    k_start=ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  stride1+ begin2+stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                            interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                            interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + stride1+begin2+stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                            +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  stride1+begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                            +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                    }

                    //k=p-1
                    else if(k<p){
                        d = data + begin1 +  stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                    }
                    //j=m-3 or m-2
                    k_start=(1+j)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  stride1+ begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + stride1+begin2+j*stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 + stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                        +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    //j=m-1 (i=1)
                    if(m%2 ==0){
                        k_start=(m%4==0)?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  stride1+begin2+(m-1)*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 + stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                        }
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        k_start=(i+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  i*stride1+begin2+j*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin1 +i* stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    k_start=(i+1)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  i*stride1+ begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + i*stride1+begin2+stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                        +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  i*stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                        +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 +  i*stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                    }
                    //j=m-3 or m-2
                    k_start=(i+j)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  i*stride1+ begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + i*stride1+begin2+j*stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                        +coeff_z*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  i*stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 + i*stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        k_start=(i+m-1)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  i*stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  i*stride1+begin2+(m-1)*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_z_xz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 + i*stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_z_xz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + i*stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_2(*(d - stride3x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(n-1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + (n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                                +coeff_z_yz*interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 +(n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x2), *(d - stride2), *(d + stride2), *(d + stride3x2)),mode);
                            }
                            
                        }
                        //j=1
                        k_start=(n%4==0)?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + (n-1)*stride1+ begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + (n-1)*stride1+begin2+stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2))  
                                                                            +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 +  (n-1)*stride1+begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) 
                                                                            +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + (n-1)*stride1+begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        k_start=(n-1+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + (n-1)*stride1+ begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                            +coeff_z_yz*interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)) 
                                                                            +coeff_z_yz*interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,interp_quad_2(*(d - stride3x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(n+m-2)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 +  (n-1)*stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 +  (n-1)*stride1+begin2+(m-1)*stride2+begin3+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin1 + (n-1)*stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride3x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin1 + (n-1)*stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }


                    //}
                    
                    //second half (adj)
                    ks1=5;
                    ks2=3;

                    for (i = 3; i + 3 < n; i += 2) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y*interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) 
                                                                                +coeff_z*interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ,mode);
                            }
                        
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_y_xy*interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ,mode);                            }

                        }
                        //j=1
                        k_start=(i+1)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + i* stride1+ begin2+stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                            +coeff_z_xz*interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ,mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + i* stride1+begin2+stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                        }

                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin1 + i* stride1+begin2+stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                        }

                        //j=m-3 or m-2 or m-1
                        while(j<m){
                            k_start=(i+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + i* stride1+ begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1))
                                                                                +coeff_z_xz*interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin1 + i* stride1+begin2+j*stride2+begin3+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1), *(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                            j+=2;
                        }
                    }
                    //i=1
                    for(j=3;j+3<m;j+=2){
                        k_start=(1+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  stride1+begin2+j*stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  stride1+begin2+j*stride2+begin3+stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if (k<p){
                            d = data + begin1 + stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
    
                    }
                    //j=1 (i+j=2)
                    k_start=ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  stride1+ begin2+stride2+begin3+k*stride3;
                        /*
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                            interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                            interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        */
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + stride1+begin2+stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                            +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                            +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);//bug when i m or p<=4, all the following quads has this problem
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  stride1+begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                            +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                            +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                    }

                    //k=p-1
                    else if(k<p){
                        d = data + begin1 +  stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y_xy*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                    }
                    //j=m-3 or m-2
                    k_start=(1+j)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  stride1+ begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + stride1+begin2+j*stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 + stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                        +coeff_y_xy*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    //j=m-1 (i=1)
                    if(m%2 ==0){
                        k_start=(m%4==0)?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                            /*
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)), 
                                interp_quad_3(*(d - stride5x2), *(d - stride3x2), *(d - stride2)),
                                interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                            */
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  stride1+begin2+(m-1)*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                            +coeff_z_xz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 + stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1))
                                                                                +coeff_z_xz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),mode);
                        }
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        k_start=(i+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2))
                                                                            +coeff_z_yz*interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  i*stride1+begin2+j*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin1 +i* stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    k_start=(i+1)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  i*stride1+ begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + i*stride1+begin2+stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                        +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  i*stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                        +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 +  i*stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)),mode);
                    }
                    //j=m-3 or m-2
                    k_start=(i+j)%4==0?ks1:ks2;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  i*stride1+ begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + i*stride1+begin2+j*stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                        +coeff_z*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  i*stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))  
                                                                        +coeff_z*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 + i*stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, coeff_x_xy*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                        +coeff_y_xy*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        k_start=(i+m-1)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  i*stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  i*stride1+begin2+(m-1)*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_z_xz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 + i*stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_x_xz*interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1))
                                                                            +coeff_z_xz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + i*stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(n-1+j)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + (n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) 
                                                                                +coeff_z_yz*interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 +(n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2), *(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            
                        }
                        //j=1
                        k_start=(n%4==0)?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + (n-1)*stride1+ begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + (n-1)*stride1+begin2+stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2))  
                                                                            +coeff_z_yz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ,mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 +  (n-1)*stride1+begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) 
                                                                            +coeff_z_yz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + (n-1)*stride1+begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        k_start=(n-1+j)%4==0?ks1:ks2;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + (n-1)*stride1+ begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                            +coeff_z_yz*interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, coeff_y_yz*interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) 
                                                                            +coeff_z_yz*interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ,mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(n+m-2)%4==0?ks1:ks2;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 +  (n-1)*stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_cubic_adj(*(d - stride3x3), *(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 +  (n-1)*stride1+begin2+(m-1)*stride2+begin3+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin1 + (n-1)*stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin1 + (n-1)*stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d,  lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }

                    /*

                    for (i = 3; i + 3 < n; i += 2) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(i+j)%4==0?5:3;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),
                                    interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) , 
                                    interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),
                                    interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ),mode);
                            }
                        
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),
                                    interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ),mode);
                            }

                        }
                        //j=1
                        k_start=(i+1)%4==0?5:3;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + i* stride1+ begin2+stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),
                                interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + i* stride1+begin2+stride2+begin3+stride3;
                           
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                        }

                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin1 + i* stride1+begin2+stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                        }

                        //j=m-3 or m-2 or m-1
                        while(j<m){
                            k_start=(i+j)%4==0?5:3;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + i* stride1+ begin2+j*stride2+begin3+k*stride3;
                               
                                predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),
                                interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ),mode);
                            }
                            //k=1
                            if(k_start==5){
                                 d = data + begin1 + i* stride1+begin2+j*stride2+begin3+stride3;
                                
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                           
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x1),*(d - stride2x1), *(d - stride1), *(d + stride1), *(d + stride2x1), *(d + stride3x1)),mode);
                            }
                            j+=2;
                        }
                    }
                    //i=1
                    for(j=3;j+3<m;j+=2){
                        k_start=(1+j)%4==0?5:3;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  stride1+begin2+j*stride2+begin3+k*stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) , 
                                    interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ) ,mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  stride1+begin2+j*stride2+begin3+stride3;
                            
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ,mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if (k<p){
                            d = data + begin1 + stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)) ,mode);
                        }
    
                    }
                    //j=1 (i+j=2)
                    for(k=3;k+3<p;k+=4){
                        d = data + begin1 +  stride1+ begin2+stride2+begin3+k*stride3;
                        
                        //predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                        //    interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)),
                        //    interp_cubic(*(d - stride3x3), *(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);
                        
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    
                    //k=1
                    //d = data + begin1 + stride1+begin2+stride2+begin3+stride3;
                        //predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1(*(d - stride1), *(d + stride1), *(d + stride3x1)),
                            //interp_quad_1(*(d - stride2), *(d + stride2), *(d + stride3x2)) , 
                            //interp_quad_1(*(d - stride3), *(d + stride3), *(d + stride3x3)) ),mode);//bug when i m or p<=4, all the following quads has this problem
                    
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  stride1+begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),
                                interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) , 
                                interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ),mode);
                    }

                    //k=p-1
                    else if(k<p){
                        d = data + begin1 +  stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),
                            interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ),mode);
                    }
                    //j=m-3 or m-2
                    k_start=(1+j)%4==0?5:3;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  stride1+ begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d,  interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + stride1+begin2+j*stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),
                                interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ),mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),
                                interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ),mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 + stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),
                            interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) ),mode);
                    }
                    //j=m-1 (i=1)
                    if(m%2 ==0){
                        k_start=(m%4==0)?5:3;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  stride1+begin2+(m-1)*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),
                                    interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ),mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 + stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),
                                                                    interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ),mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_1_adj(*(d - stride1), *(d + stride1), *(d + stride2x1)),mode);
                        }
                    }
                    //i= n-3 or n-2
                    for(j=3;j+3<m;j+=2){
                        k_start=(i+j)%4==0?5:3;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + i* stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)), 
                                interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + i*stride1+begin2+j*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
                        //k=p-3 or p-2 or p-1
                        if(k<p){
                            d = data + begin1 +i* stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                        }
                    }
                    //j=1
                    k_start=(i+1)%4==0?5:3;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  i*stride1+ begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + i*stride1+begin2+stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3(  interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),
                                interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) , 
                                interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ),mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  i*stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3(  interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),
                                interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) , 
                                interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ),mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 +  i*stride1+begin2+stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(  interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),
                            interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ),mode);
                    }
                    //j=m-3 or m-2
                    k_start=(i+j)%4==0?5:3;
                    for(k=k_start;k+3<p;k+=4){
                        d = data + begin1 +  i*stride1+ begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                    }
                    //k=1
                    if(k_start==5){
                        d = data + begin1 + i*stride1+begin2+j*stride2+begin3+stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),
                                interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ),mode);
                    }
                    //k=p-3 or p-2
                    if(k<p-1){
                        d = data + begin1 +  i*stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_ave3( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),
                                interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) , 
                                interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ),mode);
                    }
                    //k=p-1
                    else if(k<p){
                        d = data + begin1 + i*stride1+begin2+j*stride2+begin3+k*stride3;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),
                            interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2))),mode);
                    }
                    if(m%2 ==0){//j=m-1
                        k_start=(i+m-1)%4==0?5:3;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 +  i*stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 +  i*stride1+begin2+(m-1)*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),
                                            interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ),mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 + i*stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),
                                    interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ),mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + i*stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_quad_2_adj(*(d - stride2x1), *(d - stride1), *(d + stride1)),mode);
                        }
                    }
                    //i=n-1 (odd)
                    if (n % 2 == 0) {
                        for(j=3;j+3<m;j+=2){
                            k_start=(n-1+j)%4==0?5:3;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 + (n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)), 
                                    interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)) ),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            //k=p-3 or p-2 or p-1
                            if(k<p){
                                d = data + begin1 +(n-1)* stride1+begin2+j*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x2),*(d - stride2x2), *(d - stride2), *(d + stride2), *(d + stride2x2), *(d + stride3x2)),mode);
                            }
                            
                        }
                        //j=1
                        k_start=(n%4==0)?5:3;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + (n-1)*stride1+ begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,  interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + (n-1)*stride1+begin2+stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) , 
                                    interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ),mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 +  (n-1)*stride1+begin2+stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) , 
                                    interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ),mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + (n-1)*stride1+begin2+stride2+begin3+k*stride3;

                            predict_error+=quantize_integrated(d - data, *d, interp_quad_1_adj(*(d - stride2), *(d + stride2), *(d + stride2x2)) ,mode);
                        }
                        //j=m-3 or m-2
                        k_start=(n-1+j)%4==0?5:3;
                        for(k=k_start;k+3<p;k+=4){
                            d = data + begin1 + (n-1)*stride1+ begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                        }
                        //k=1
                        if(k_start==5){
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2+begin3+stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) , 
                                    interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)) ),mode);
                        }
                        //k=p-3 or p-2
                        if(k<p-1){
                            d = data + begin1 +  (n-1)*stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)) , 
                                    interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)) ),mode);
                        }
                        //k=p-1
                        else if(k<p){
                            d = data + begin1 + (n-1)*stride1+begin2+j*stride2+begin3+k*stride3;
                            predict_error+=quantize_integrated(d - data, *d,interp_quad_2_adj(*(d - stride2x2), *(d - stride2), *(d + stride2)),mode);
                        }
                        if(m%2 ==0){//j=m-1
                            k_start=(n+m-2)%4==0?5:3;
                            for(k=k_start;k+3<p;k+=4){
                                d = data + begin1 +  (n-1)*stride1+ begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_cubic_adj(*(d - stride3x3),*(d - stride2x3), *(d - stride3), *(d + stride3), *(d + stride2x3), *(d + stride3x3)),mode);
                            }
                            //k=1
                            if(k_start==5){
                                d = data + begin1 +  (n-1)*stride1+begin2+(m-1)*stride2+begin3+stride3;
                                predict_error+=quantize_integrated(d - data, *d,interp_quad_1_adj(*(d - stride3), *(d + stride3), *(d + stride2x3)),mode);
                            }
                            //k=p-3 or p-2
                            if(k<p-1){
                                d = data + begin1 + (n-1)*stride1+begin2+(m-1)*stride2+begin3+k*stride3;
                                predict_error+=quantize_integrated(d - data, *d, interp_quad_2_adj(*(d - stride2x3), *(d - stride3), *(d + stride3)),mode);
                            }
                            //k=p-1
                            else if(k<p){
                                d = data + begin1 + (n-1)*stride1+begin2+(m-1)*stride2+begin3+(p-1)*stride3;
                                predict_error+=quantize_integrated(d - data, *d, lorenzo_3d(*(d-stride1-stride2-stride3),*(d-stride1-stride2),*(d-stride1-stride3),*(d-stride1),*(d-stride2-stride3),*(d-stride2),*(d-stride3)),mode);
                            }
                        }
                    }
                    */
                }
            }
            return predict_error;
        }

        double block_interpolation_2d_cross(T *data, size_t begin1, size_t end1, size_t begin2, size_t end2, size_t stride1,size_t stride2, const std::string &interp_func, const PredictorBehavior pb,const QoZ::Interp_Meta &meta,int tuning=0) {
           // std::cout<<"cst"<<std::endl;
            size_t n = (end1 - begin1) / stride1 + 1;
            if (n <= 1) {
                return 0;
            }
            size_t m = (end2 - begin2) / stride2 + 1;
            if (m <= 1) {
                return 0;
            }

            double predict_error = 0;
            size_t stride3x1=3*stride1,stride3x2=3*stride2,stride5x1=5*stride1,stride5x2=5*stride2;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            if (interp_func == "linear"|| n<5 || m<5 ) {//nmcond temp added
                
                for (size_t i = 1; i + 1 < n; i += 2) {
                    for(size_t j=1;j+1<m;j+=2){
                        T *d = data + begin1 + i* stride1+begin2+j*stride2;
                       
                        predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);

                    }
                    if(m%2 ==0){
                        T *d = data + begin1 + i * stride1+begin2+(m-1)*stride2;

                        if(i<3 or i+3>=n or m<4)
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1 - stride2), *(d + stride1 - stride2)),mode);//this is important. Not sure whether it is good.
                        else
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_linear1(*(d - stride3x1 - stride3x2),*(d - stride1 - stride2))
                                                            , interp_linear1(*(d + stride3x1 - stride3x2),*(d + stride1 - stride2))),mode);//this is important. Not sure whether it is good.
                    }
                }
                if (n % 2 == 0) {
                    for(size_t j=1;j+1<m;j+=2){

                        T *d = data + begin1 + (n-1) * stride1+begin2+j*stride2;
                        if(n<4 or j<3 or j+3>=m)
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1 - stride2), *(d - stride1 + stride2)),mode);//this is important. Not sure whether it is good.
                        else
                            predict_error+=quantize_integrated(d - data, *d, interp_linear(interp_linear1(*(d - stride3x1 - stride3x2),*(d - stride1 - stride2))
                                                            , interp_linear1(*(d - stride3x1 + stride3x2),*(d - stride1 + stride2))),mode);//this is important. Not sure whether it is good.
                    }
                    if(m%2 ==0){
                        T *d = data + begin1 + (n-1) * stride1+begin2+(m-1)*stride2;
                        if(n<4 or m<4)
                            predict_error+=quantize_integrated(d - data, *d, *(d - stride1 - stride2),mode);//this is important. Not sure whether it is good.
                        else
                            predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d - stride3x1 - stride3x2),*(d - stride1 - stride2) ),mode);//this is important. Not sure whether it is good.
                    }          
                }
                    
            }
            else{//cubic
                //adaptive todo
                auto interp_cubic=meta.cubicSplineType==0?interp_cubic_1<T>:interp_cubic_2<T>;
                size_t i,j;
                T *d;
                for (i = 3; i + 3 < n; i += 2) {
                    for(j=3;j+3<m;j+=2){
                        d = data + begin1 + i* stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1-stride3x2), *(d - stride1-stride2), *(d + stride1+stride2), *(d + stride3x1+stride3x2))
                                                        ,interp_cubic(*(d +stride3x1- stride3x2), *(d +stride1- stride2), *(d -stride1+ stride2), *(d -stride3x1+ stride3x2)) ),mode);
                    }
                    //j=1
                    d = data + begin1 + i* stride1+ begin2+stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1( *(d - stride1-stride2),*(d + stride1+stride2),*(d + stride3x1+stride3x2) )
                                                    ,interp_quad_1( *(d + stride1-stride2),*(d - stride1+stride2),*(d - stride3x1+stride3x2) ) ),mode);
                                       
                    //j=m-3 or m-2
                    d = data +begin1 + i* stride1+ begin2+j*stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2( *(d - stride3x1-stride3x2),*(d - stride1-stride2),*(d + stride1+stride2) )
                                                    ,interp_quad_2( *(d + stride3x1-stride3x2),*(d + stride1-stride2),*(d - stride1+stride2) ) ),mode);
                    
                    //j=m-1
                    if(m%2 ==0){
                        d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                        if(i>=5 and i+5<n)
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_3( *(d - stride5x1-stride5x2),*(d - stride3x1-stride3x2),*(d - stride1-stride2) )
                                                   ,interp_quad_2( *(d + stride5x1-stride5x2),*(d + stride3x1-stride3x2),*(d + stride1-stride2) ) ),mode);
                        else
                            predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_linear1( *(d - stride3x1-stride3x2),*(d - stride1-stride2) )
                                                    ,interp_linear1(*(d + stride3x1-stride3x2),*(d + stride1-stride2) ) ),mode);
                    }
                }
               
                //i=1
                for(j=3;j+3<m;j+=2){
                    d = data + begin1 + stride1+begin2+j*stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_1( *(d - stride1-stride2),*(d + stride1+stride2),*(d + stride3x1+stride3x2) )
                                                    ,interp_quad_1( *(d - stride1+stride2),*(d + stride1-stride2),*(d + stride3x1-stride3x2) ) ),mode);
                }
                //j=1

                d = data + begin1 + stride1+ begin2+stride2;
                //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);//2d linear
                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad1( *(d - stride1-stride2),*(d + stride1+stride2),*(d + stride3x1+stride3x2) )
                //                                    ,interp_linear( *(d + stride1-stride2),*(d - stride1+stride2) ) ),mode);//2d linear+quad
                predict_error+=quantize_integrated(d - data, *d, interp_quad_1( *(d - stride1-stride2),*(d + stride1+stride2),*(d + stride3x1+stride3x2) ),mode);//1d quad

                //j=m-3 or m-2
                d = data +begin1 + stride1+ begin2+j*stride2;

                //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);//2d linear
                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad2( *(d + stride3x1-stride3x2),*(d + stride1-stride2),*(d - stride1+stride2) )
                //                                    ,interp_linear( *(d - stride1-stride2),*(d + stride1+stride2) ) ),mode);//2d linear+quad
                predict_error+=quantize_integrated(d - data, *d, interp_quad_2( *(d + stride3x1-stride3x2),*(d + stride1-stride2),*(d - stride1+stride2) ),mode);//1d quad
                
                //j=m-1
                if(m%2 ==0){
                    d = data + begin1 + stride1+begin2+(m-1)*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1-stride2), *(d + stride1-stride2)),mode);//1d linear
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d + stride3x1-stride3x2), *(d + stride1-stride2)),mode);//1d cross linear
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1-stride2), *(d + stride1-stride2),*(d + stride3x1-stride2)),mode);//1d quad

                }

                //i= n-3 or n-2
                for(j=3;j+3<m;j+=2){
                   
                    d = data + begin1 + i*stride1+begin2+j*stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad_2( *(d - stride3x1-stride3x2),*(d - stride1-stride2),*(d + stride1+stride2) )
                                                    ,interp_quad_2( *(d - stride3x1+stride3x2),*(d - stride1+stride2),*(d + stride1-stride2) ) ),mode);

                }
                //j=1
                d = data + begin1 + i*stride1+ begin2+stride2;
                //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);//2d linear
                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad1( *(d + stride1-stride2),*(d - stride1+stride2),*(d - stride3x1+stride3x2) )
                //                                    ,interp_linear( *(d - stride1-stride2),*(d + stride1+stride2) ) ),mode);//2d linear+quad
                predict_error+=quantize_integrated(d - data, *d, interp_quad_1( *(d + stride1-stride2),*(d - stride1+stride2),*(d - stride3x1+stride3x2) ),mode);//1d quad
                
                //j=m-3 or m-2
                d = data +begin1 + i*stride1+ begin2+j*stride2;
           
                //predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 - stride2), *(d + stride1 + stride2),*(d - stride1 + stride2), *(d + stride1 - stride2)),mode);//2d linear
                //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_quad2( *(d - stride3x1-stride3x2),*(d - stride1-stride2),*(d + stride1+stride2) )
                //                                    ,interp_linear( *(d + stride1-stride2),*(d - stride1+stride2) ) ),mode);//2d linear+quad
                predict_error+=quantize_integrated(d - data, *d, interp_quad_2( *(d - stride3x1-stride3x2),*(d - stride1-stride2),*(d + stride1+stride2) ),mode);//1d quad
                
                //j=m-1
                if(m%2 ==0){
                    d = data + begin1 + i * stride1+begin2+(m-1)*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1-stride2), *(d + stride1-stride2)),mode);//1d linear
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear1(*(d + stride3x1-stride3x2), *(d + stride1-stride2)),mode);//1d cross linear
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d + stride1-stride2), *(d - stride1-stride2),*(d - stride3x1-stride2)),mode);//1d quad
                }

                //i=n-1 (odd)
                if (n % 2 == 0) {
                    for(j=3;j+3<m;j+=2){
                        d = data + begin1 + (n-1)*stride1+begin2+j*stride2;
                        //predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_linear1(*(d - stride3x1-stride3x2), *(d - stride1-stride2)) ,interp_linear1(*(d - stride3x1+stride3x2), *(d - stride1+stride2)) ),mode);//2d cross
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride1-stride3x2), *(d -  stride1-stride2), *(d - stride1+ stride2), *(d - stride1+ stride3x2)),mode);//1d cubic


                    }
                    //j=1
                    d = data + begin1 + (n-1)*stride1+ begin2+stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear1( *(d - stride3x1+stride3x2), *(d - stride1+stride2)),mode);//1d linear
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_1(*(d - stride1-stride2), *(d - stride1+stride2),*(d - stride1+stride3x2)),mode);//1d quad

                    //j=m-3 or m-2
                    d = data +begin1 + (n-1)*stride1+ begin2+j*stride2;
                    //predict_error+=quantize_integrated(d - data, *d, interp_linear1( *(d - stride3x1-stride3x2), *(d - stride1-stride2)),mode);//1d linear
                    predict_error+=quantize_integrated(d - data, *d, interp_quad_2(*(d - stride1-stride3x2), *(d - stride1-stride2),*(d - stride1+stride2)),mode);//1d quad
                    //j=m-1
                    if(m%2 ==0){
                        d = data + begin1 + (n-1) * stride1+begin2+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_quad_3(*(d - stride5x1-stride5x2), *(d - stride3x1-stride3x2), *(d - stride1-stride2)),mode);
                    }
                }         
            }
            return predict_error;
        }

        double block_interpolation_2d_aftercross(T *data, size_t begin1, size_t end1, size_t begin2, size_t end2, size_t stride1,size_t stride2, const std::string &interp_func, const PredictorBehavior pb,const QoZ::Interp_Meta &meta,int tuning=0) {
            size_t n = (end1 - begin1) / stride1 + 1;
            
            size_t m = (end2 - begin2) / stride2 + 1;
            if (n<=1&& m <= 1) {
                return 0;
            }
            double predict_error = 0;
            size_t stride3x1=3*stride1,stride3x2=3*stride2,stride5x1=5*stride1,stride5x2=5*stride2;
            int mode=(pb == PB_predict_overwrite)?tuning:-1;
            if (interp_func == "linear"|| n<5 || m<5 ) {//nmcond temp added
                size_t i,j;
                for (i = 1; i + 1 < n; i += 1) {
                    for(j=1+(i%2);j+1<m;j+=2){
                        T *d = data + begin1 + i* stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_2d(*(d - stride1 ), *(d + stride1 ),*(d  + stride2), *(d  - stride2)),mode);

                    }

                    //j=0
                    if(i%2==1 and begin2==0){
                        T *d = data + begin1 + i* stride1+begin2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1 ), *(d + stride1 )),mode);
                    }
                    //j=m-1, j wont be 0
                    if(j==m-1){
                        T *d = data + begin1 + i* stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d - stride1 ), *(d + stride1 )),mode);
                    }
                }
                //i=0
                if(begin1==0){
                    for(j=1;j+1<m;j+=2){
                        T *d = data + begin1 +begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d  + stride2), *(d  - stride2)),mode);
                    }
                
                //j=m-1, j wont be 0
                    if(j==m-1){
                        T *d = data + begin1 +begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, *(d-stride2),mode);//for simplicity,may extend to 2d.
                    }
                }
                //i=n-1
                if(n>1){
                    for(j=1+(n-1)%2;j+1<m;j+=2){
                        T *d = data + begin1 +(n-1)*stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear(*(d  + stride2), *(d  - stride2)),mode);
                    }
                    //j=0
                    if((n-1)%2==1 and begin2==0){
                        
                        T *d = data + begin1 +(n-1)*stride1+begin2;

                        predict_error+=quantize_integrated(d - data, *d, *(d-stride1),mode);//for simplicity,may extend to 2d.
                    }
                    //j=m-1, j wont be 0
                    if( j==m-1){
                        
                        T *d = data + begin1 +(n-1)*stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, *(d-stride1),mode);//for simplicity,may extend to 2d.
                    }
                }
                    
            }
            else{//cubic
                //adaptive todo
                auto interp_cubic=meta.cubicSplineType==0?interp_cubic_1<T>:interp_cubic_2<T>;
                size_t i,j;
                T *d;
                for (i = 3; i + 3 < n; i += 1) {
                    for(j=3+(i%2);j+3<m;j+=2){
                        d = data + begin1 + i* stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d, interp_linear( interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1))
                                                        ,interp_cubic(*(d - stride3x2), *(d - stride2), *(d+ stride2), *(d + stride3x2)) ),mode);
                    }
                    //j=0
                    if(i%2==1 and begin2==0){
                        d = data + begin1 + i* stride1+begin2;
                        predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    }
                    //j=1 or 2 
                    d = data + begin1 + i* stride1+begin2+(1+(i%2))*stride2;
                    predict_error+=quantize_integrated(d - data, *d, interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    
                    //j=m-3 or m-2, j wont be 2.
                    d = data + begin1 + i* stride1+begin2+j*stride2;
                    predict_error+=quantize_integrated(d - data, *d,  interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                   
                    //j=m-1
                    if(j+2==m-1){
                        d = data + begin1 + i* stride1+begin2+(m-1)*stride2;
                        predict_error+=quantize_integrated(d - data, *d,  interp_cubic(*(d - stride3x1), *(d - stride1), *(d + stride1), *(d + stride3x1)),mode);
                    }
                   
                }
                std::vector<size_t> boundary_is=(n>5)?std::vector<size_t>{0,1,2,n-3,n-2,n-1}:std::vector<size_t>{0,1,2,n-2,n-1};
                
                for(auto ii:boundary_is){
                    if(ii==0 and begin1!=0)
                        continue;
                    for(j=3+(ii%2);j+3<m;j+=2){
                        d = data + begin1+ii*stride1+begin2+j*stride2;
                        predict_error+=quantize_integrated(d - data, *d,interp_cubic(*(d - stride3x2), *(d - stride2), *(d+ stride2), *(d + stride3x2) ),mode);
                    }

                    std::vector<size_t> boundary_js=(ii%2)?std::vector<size_t>{0,2,j}:std::vector<size_t>{1,j};
                    if(j+2==m-1)
                        boundary_js.push_back(m-1);

                    for(auto jj:boundary_js){
                        if(begin2!=0 and jj==0)
                            continue;
                        d = data + begin1 + ii* stride1+begin2+jj*stride2;
                        T v1;
                        if(ii==0){
                            if(n>5)
                                v1=interp_quad_3(*(d + stride5x1), *(d+ stride3x1), *(d + stride1) );
                            else
                                v1=interp_linear1( *(d+ stride3x1), *(d + stride1) );
                        }
                        else if(ii==1)
                            v1=interp_quad_1(*(d - stride1), *(d+ stride1), *(d + stride3x1) );
                        else if (ii==n-2)
                            v1=interp_quad_2(*(d - stride3x1), *(d- stride1), *(d + stride1) );
                        else if (ii==n-1){
                            if(n>5)
                                v1=interp_quad_3(*(d - stride5x1), *(d- stride3x1), *(d - stride1) );
                            else
                                v1=interp_linear1( *(d- stride3x1), *(d - stride1) );
                        }
                        else{//i==2 or n-3
                            if(n==5)
                                v1=interp_linear(*(d - stride1), *(d+ stride1));
                            else if (ii==2)
                                v1=interp_quad_1(*(d - stride1), *(d+ stride1), *(d + stride3x1) );
                            else
                                v1=interp_quad_2(*(d - stride3x1), *(d- stride1), *(d + stride1) );
                        }

                        T v2;
                        if(jj==0){
                            if(m>5)
                                v2=interp_quad_3( *(d + stride5x2), *(d+ stride3x2), *(d + stride2) );
                            else
                                v2=interp_linear1( *(d+ stride3x2), *(d + stride2) );
                        }
                        else if (jj==1 or jj==2)
                            v2=interp_quad_1( *(d - stride2), *(d+ stride2), *(d + stride3x2) );
                        else if(jj==m-1){
                            if(m>5)
                                v2=interp_quad_3( *(d - stride5x2), *(d- stride3x2), *(d - stride2) );
                            else
                                v2=interp_linear1( *(d- stride3x2), *(d - stride2) );
                        }
                        else
                            v2=interp_quad_2( *(d - stride3x2), *(d- stride2), *(d + stride2) );

                        predict_error+=quantize_integrated(d - data, *d,interp_linear( v1,v2 ),mode);

                    }

                }
            }
            return predict_error;
        }
        
        template<uint NN = N>
        typename std::enable_if<NN == 1, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const QoZ::Interp_Meta & meta, size_t stride = 1,int tuning=0,size_t cross_block=0,int regressive=0) {//regressive to reduce into meta.
            return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb,meta,tuning);
        }


        template<uint NN = N>
        typename std::enable_if<NN == 2, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const QoZ::Interp_Meta & meta, size_t stride = 1,int tuning=0,size_t cross_block=0,int regressive=0) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            //bool full_adjacent_interp=false;
            uint8_t paradigm=meta.interpParadigm;
            uint8_t direction=meta.interpDirection;
            assert(direction<2);
            if(paradigm==0){
                const std::array<int, N> dims = dimension_sequences[direction];
                
                
                if(!regressive or stride!=1){
                    
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                                stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning);
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                    }
                }
                else{
                    std::vector<size_t> block_begin(begin.begin(),begin.end());
                    std::vector<size_t> block_end(end.begin(),end.end());
                    std::vector<size_t> sparsity={2,2};
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        std::vector<size_t> cur_begin=block_begin;
                        cur_begin[dims[1]]=j;
                        std::vector<size_t> cur_end=block_end;
                        cur_end[dims[1]]=j;
                        predict_error += block_interpolation_1d_regressive(data,block_begin,block_end,dims[0],cur_begin,cur_end,sparsity,
                                                        stride,interp_func,pb,meta,tuning);
                    }
                    sparsity[dims[0]]=1;
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                        std::vector<size_t> cur_begin=block_begin;
                        cur_begin[dims[0]]=i;
                        std::vector<size_t> cur_end=block_end;
                        cur_end[dims[0]]=i;
                        predict_error += block_interpolation_1d_regressive(data,block_begin,block_end,dims[1],cur_begin,cur_end,sparsity,
                                                        stride,interp_func,pb,meta,tuning);
                    }

                }
            }
            
            else if(paradigm<3){//md or hd
                const std::array<int, N> dims = dimension_sequences[0];
                std::array<double,3>dim_coeffs=meta.dimCoeffs;
                
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                    size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning);
                }
                for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride2x : 0); i <= end[dims[0]]; i += stride2x) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                }
                size_t begin_offset1=begin[dims[0]]*dimension_offsets[dims[0]];
                size_t begin_offset2=begin[dims[1]]*dimension_offsets[dims[1]];
                //size_t stride1=dimension_offsets[dims[0]]
                predict_error+=block_interpolation_2d(data, begin_offset1,
                                                            begin_offset1 +
                                                            (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                            begin_offset2,
                                                            begin_offset2 +
                                                            (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb,std::array<double,2>{dim_coeffs[0],dim_coeffs[1]},meta,tuning);
            }
           
            else if(paradigm==3){//cross
                const std::array<int, N> dims = dimension_sequences[0];
                size_t begin_offset1=begin[dims[0]]*dimension_offsets[dims[0]];
                size_t begin_offset2=begin[dims[1]]*dimension_offsets[dims[1]];
                predict_error+=block_interpolation_2d_cross(data, begin_offset1,
                                                            begin_offset1 +
                                                            (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                            begin_offset2,
                                                            begin_offset2 +
                                                            (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                predict_error+=block_interpolation_2d_aftercross(data, begin_offset1,
                                                            begin_offset1 +
                                                            (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                            begin_offset2,
                                                            begin_offset2 +
                                                            (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
            }

            return predict_error;
        }




        template<uint NN = N>
        typename std::enable_if<NN == 3, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const QoZ::Interp_Meta & meta, size_t stride = 1,int tuning=0,size_t cross_block=0,int regressive=0) {//cross block: 0 or conf.num

            double predict_error = 0;
            size_t stride2x = stride * 2;
           //bool full_adjacent_interp=false;
            uint8_t paradigm=meta.interpParadigm;
            uint8_t direction=meta.interpDirection;
            bool fallback_2d=direction>=6;
            if (fallback_2d){
                direction-=6;
            }
            assert(direction<6);
            if(paradigm==0){
                const std::array<int, N> dims = dimension_sequences[direction];
                //if (cross_block==0){
                if(!fallback_2d){
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[0]] - begin[dims[0]]) *
                                                                    dimension_offsets[dims[0]],
                                                                    stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  begin[dims[2]] * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning);
                        }
                    }
                }
                else{
                    //dims[0] fused.
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + 1 : 0); i <= end[dims[0]]; i += 1) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + 1 : 0); i <= end[dims[0]]; i += 1) {
                        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  begin[dims[2]] * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning);
                        }
                    }
                }
                //}
                    /*
                else{
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                           // std::cout<<"phase 1"<<" "<<j<<" "<<k<<std::endl;
                            predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[0]] - begin[dims[0]]) *
                                                                    dimension_offsets[dims[0]],
                                                                    stride * dimension_offsets[dims[0]], interp_func, pb,tuning,2,begin[dims[0]],stride,dims[0]);
                        }
                    }
                   // size_t iidx=begin[dims[0]] ? 1: 0;
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            if(i%(stride2x)==0){
                                //if(tuning==0 and stride==1)
                                   // std::cout<<"phase 2"<<" "<<i<<" "<<k<<std::endl;
                                predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,tuning,2,begin[dims[1]],stride,dims[1]);
                            }
                            else{
                                //std::cout<<"phase 3"<<" "<<i<<" "<<k<<std::endl;
                                predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,tuning,1,begin[dims[1]],stride,dims[1]);
                            }
                           
                        }
                        //iidx++;
                    }
                    //iidx=begin[dims[0]] ? 1: 0;
                    //size_t jjdx=begin[dims[1]] ? 1: 0;
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                         for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  begin[dims[2]] * dimension_offsets[dims[2]];
                            if(i%(stride2x)==0 and j%(stride2x)==0){
                                //if(tuning==0 and stride==1)
                                    //std::cout<<"phase 4"<<" "<<i<<" "<<j<<std::endl;
                                predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                        begin_offset +
                                                                        (end[dims[2]] - begin[dims[2]]) *
                                                                        dimension_offsets[dims[2]],
                                                                        stride * dimension_offsets[dims[2]], interp_func, pb,tuning,2,begin[dims[2]],stride,dims[2]);
                            }
                            else{
                                //std::cout<<"phase 5"<<" "<<i<<" "<<j<<std::endl;
                                predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                        begin_offset +
                                                                        (end[dims[2]] - begin[dims[2]]) *
                                                                        dimension_offsets[dims[2]],
                                                                        stride * dimension_offsets[dims[2]], interp_func, pb,tuning,1,begin[dims[2]],stride,dims[2]);
                            }


                            
                            //jjdx++;

                        }
                        //iidx++;
                    }

                }
                */
            }
            
            else if (paradigm==1){
                if(!fallback_2d){
                    const std::array<int, N> dims = dimension_sequences[0];
                    std::array<double,3>dim_coeffs=meta.dimCoeffs;
                    //std::cout<<dim_coeffs[0]<<std::endl;
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[0]] - begin[dims[0]]) *
                                                                    dimension_offsets[dims[0]],
                                                                    stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride2x : 0); i <= end[dims[0]]; i += stride2x) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride2x : 0); i <= end[dims[0]]; i += stride2x) {
                        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  begin[dims[2]] * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning);
                        }
                    }
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset1 = begin[dims[0]] * dimension_offsets[dims[0]] + k * dimension_offsets[dims[2]];
                            size_t begin_offset2 =  begin[dims[1]] * dimension_offsets[dims[1]];
                                                
                            predict_error += block_interpolation_2d(data, begin_offset1,
                                                                    begin_offset1 +
                                                                    (end[dims[0]] - begin[dims[0]]) *
                                                                    dimension_offsets[dims[0]],
                                                                    begin_offset2,
                                                                    begin_offset2 +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[0]], stride * dimension_offsets[dims[1]],interp_func, pb,std::array<double,2>{dim_coeffs[dims[0]],dim_coeffs[dims[1]]},meta,tuning);
                        }


                        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                            size_t begin_offset1 = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
                            size_t begin_offset2 =  begin[dims[2]] * dimension_offsets[dims[2]];
                                                
                            predict_error += block_interpolation_2d(data, begin_offset1,
                                                                    begin_offset1 +
                                                                    (end[dims[0]] - begin[dims[0]]) *
                                                                    dimension_offsets[dims[0]],
                                                                    begin_offset2,
                                                                    begin_offset2 +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[0]], stride * dimension_offsets[dims[2]],interp_func, pb,std::array<double,2>{dim_coeffs[dims[0]],dim_coeffs[dims[2]]},meta,tuning);
                        }
                        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride2x : 0); i <= end[dims[0]]; i += stride2x) {
                            size_t begin_offset1 = begin[dims[1]] * dimension_offsets[dims[1]] + i * dimension_offsets[dims[0]];
                            size_t begin_offset2 =  begin[dims[2]] * dimension_offsets[dims[2]];
                                                
                            predict_error += block_interpolation_2d(data, begin_offset1,
                                                                    begin_offset1 +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    begin_offset2,
                                                                    begin_offset2 +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[1]], stride * dimension_offsets[dims[2]],interp_func, pb,std::array<double,2>{dim_coeffs[dims[1]],dim_coeffs[dims[2]]},meta,tuning);
                        }
                        size_t begin_offset1 = begin[dims[0]] * dimension_offsets[dims[0]] ;
                        size_t begin_offset2 = begin[dims[1]] * dimension_offsets[dims[1]] ;
                        size_t begin_offset3 =  begin[dims[2]] * dimension_offsets[dims[2]];
                        predict_error += block_interpolation_3d(data, begin_offset1,
                                                                    begin_offset1 +
                                                                    (end[dims[0]] - begin[dims[0]]) *
                                                                    dimension_offsets[dims[0]],
                                                                    begin_offset2,
                                                                    begin_offset2 +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    begin_offset3,
                                                                    begin_offset3 +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[0]],stride * dimension_offsets[dims[1]], stride * dimension_offsets[dims[2]],interp_func,pb,dim_coeffs,meta,tuning);
                }
                else{
                    const std::array<int, N> dims = dimension_sequences[direction];
                    std::array<double,3>dim_coeffs=meta.dimCoeffs;
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + 1 : 0); i <= end[dims[0]]; i += 1) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + 1: 0); i <= end[dims[0]]; i += 1) {
                        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  begin[dims[2]] * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + 1 : 0); i <= end[dims[0]]; i += 1) {
                            size_t begin_offset1 = begin[dims[1]] * dimension_offsets[dims[1]] + i * dimension_offsets[dims[0]];
                            size_t begin_offset2 =  begin[dims[2]] * dimension_offsets[dims[2]];
                                                
                            predict_error += block_interpolation_2d(data, begin_offset1,
                                                                    begin_offset1 +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    begin_offset2,
                                                                    begin_offset2 +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[1]], stride * dimension_offsets[dims[2]],interp_func, pb,std::array<double,2>{dim_coeffs[dims[1]],dim_coeffs[dims[2]]},meta,tuning);
                        }

                }
            }
            else if (paradigm==2){
                const std::array<int, N> dims = dimension_sequences[direction];
                std::array<double,3>dim_coeffs=meta.dimCoeffs;
                //don't do md interp on dims[0]
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                        size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              k * dimension_offsets[dims[2]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning);
                    }
                }
                for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                              k * dimension_offsets[dims[2]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                    }
                }
                for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              begin[dims[2]] * dimension_offsets[dims[2]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning);
                    }
                }
                    

                for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                    size_t begin_offset1 = begin[dims[1]] * dimension_offsets[dims[1]] + i * dimension_offsets[dims[0]];
                    size_t begin_offset2 =  begin[dims[2]] * dimension_offsets[dims[2]];
                                            
                    predict_error += block_interpolation_2d(data, begin_offset1,
                                                            begin_offset1 +
                                                            (end[dims[1]] - begin[dims[1]]) *
                                                            dimension_offsets[dims[1]],
                                                            begin_offset2,
                                                            begin_offset2 +
                                                            (end[dims[2]] - begin[dims[2]]) *
                                                            dimension_offsets[dims[2]],
                                                            stride * dimension_offsets[dims[1]], stride * dimension_offsets[dims[2]],interp_func, pb,std::array<double,2>{dim_coeffs[dims[1]],dim_coeffs[dims[2]]},meta,tuning);
                }
                    
            }
            return predict_error;
        }


        template<uint NN = N>
        typename std::enable_if<NN == 4, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const QoZ::Interp_Meta & meta, size_t stride = 1,int tuning=0,size_t cross_block=0,int regressive=0) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            uint8_t paradigm=meta.interpParadigm;
            uint8_t direction=meta.interpDirection;
            assert(direction<24);
            //max_error = 0;
            const std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning);
                    }
                }
            }
//            printf("%.8f ", max_error);
           // max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                    }
                }
            }
//            printf("%.8f ", max_error);
            //max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              begin[dims[2]] * dimension_offsets[dims[2]] +
                                              t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning);
                    }
                }
            }

//            printf("%.8f ", max_error);
          //  max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0); k <= end[dims[2]]; k += stride) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                begin[dims[3]] * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[3]] - begin[dims[3]]) *
                                                                dimension_offsets[dims[3]],
                                                                stride * dimension_offsets[dims[3]], interp_func, pb,meta,tuning);
                    }
                }
            }
//            printf("%.8f \n", max_error);
            return predict_error;
        }


        double block_interpolation_block3d(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func,const QoZ::Interp_Meta & meta, size_t stride = 1,int tuning=0,size_t cross_block=0) {//cross block: 0 or conf.num

            double predict_error = 0;
            size_t stride2x = stride * 2;
            uint8_t paradigm=meta.interpParadigm;
            uint8_t direction=meta.interpDirection;
            assert(direction<6);
            //if(direction!=6){
            const std::array<int, N> dims = dimension_sequences[direction];

                if (cross_block==0){
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[0]] - begin[dims[0]]) *
                                                                    dimension_offsets[dims[0]],
                                                                    stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  begin[dims[2]] * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[2]] - begin[dims[2]]) *
                                                                    dimension_offsets[dims[2]],
                                                                    stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning);
                        }
                    }
                }
                    
                else{

                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[0]] - begin[dims[0]]) *
                                                                    dimension_offsets[dims[0]],
                                                                    stride * dimension_offsets[dims[0]], interp_func, pb,meta,tuning,2,begin[dims[0]],stride,dims[0]);
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                        for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                                  k * dimension_offsets[dims[2]];
                            if(i%(stride2x)==0){
                                predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning,2,begin[dims[1]],stride,dims[1]);
                            }
                            else{
                                predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                    begin_offset +
                                                                    (end[dims[1]] - begin[dims[1]]) *
                                                                    dimension_offsets[dims[1]],
                                                                    stride * dimension_offsets[dims[1]], interp_func, pb,meta,tuning,1,begin[dims[1]],stride,dims[1]);
                            }
                           
                        }
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                         for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                            size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                                  begin[dims[2]] * dimension_offsets[dims[2]];
                            if(i%(stride2x)==0 and j%(stride2x)==0){
                                predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                        begin_offset +
                                                                        (end[dims[2]] - begin[dims[2]]) *
                                                                        dimension_offsets[dims[2]],
                                                                        stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning,2,begin[dims[2]],stride,dims[2]);
                            }
                            else{ 
                                predict_error += block_interpolation_1d_cross(data, begin_offset,
                                                                        begin_offset +
                                                                        (end[dims[2]] - begin[dims[2]]) *
                                                                        dimension_offsets[dims[2]],
                                                                        stride * dimension_offsets[dims[2]], interp_func, pb,meta,tuning,1,begin[dims[2]],stride,dims[2]);
                            } 
                        }
                    }
                }
            //}
            /*
            else{
                const std::array<int, N> dims = dimension_sequences[0];
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                        size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              k * dimension_offsets[dims[2]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                stride * dimension_offsets[dims[0]], interp_func, pb,tuning);
                    }
                }
                for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride2x : 0); i <= end[dims[0]]; i += stride2x) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                              k * dimension_offsets[dims[2]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[1]], interp_func, pb,tuning);
                    }
                }
                for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride2x : 0); i <= end[dims[0]]; i += stride2x) {
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              begin[dims[2]] * dimension_offsets[dims[2]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[2]], interp_func, pb,tuning);
                    }
                }
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                        size_t begin_offset1 = begin[dims[0]] * dimension_offsets[dims[0]] + k * dimension_offsets[dims[2]];
                        size_t begin_offset2 =  begin[dims[1]] * dimension_offsets[dims[1]];
                                            
                        predict_error += block_interpolation_2d(data, begin_offset1,
                                                                begin_offset1 +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                begin_offset2,
                                                                begin_offset2 +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[0]], stride * dimension_offsets[dims[1]],interp_func, pb,tuning);
                    }
                    for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                        size_t begin_offset1 = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
                        size_t begin_offset2 =  begin[dims[2]] * dimension_offsets[dims[2]];
                                            
                        predict_error += block_interpolation_2d(data, begin_offset1,
                                                                begin_offset1 +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                begin_offset2,
                                                                begin_offset2 +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[0]], stride * dimension_offsets[dims[2]],interp_func, pb,tuning);
                    }
                    for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride2x : 0); i <= end[dims[0]]; i += stride2x) {
                        size_t begin_offset1 = begin[dims[1]] * dimension_offsets[dims[1]] + i * dimension_offsets[dims[0]];
                        size_t begin_offset2 =  begin[dims[2]] * dimension_offsets[dims[2]];
                                            
                        predict_error += block_interpolation_2d(data, begin_offset1,
                                                                begin_offset1 +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                begin_offset2,
                                                                begin_offset2 +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[1]], stride * dimension_offsets[dims[2]],interp_func, pb,tuning);
                    }
                    size_t begin_offset1 = begin[dims[0]] * dimension_offsets[dims[0]] ;
                    size_t begin_offset2 = begin[dims[1]] * dimension_offsets[dims[1]] ;
                    size_t begin_offset3 =  begin[dims[2]] * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_3d(data, begin_offset1,
                                                                begin_offset1 +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                begin_offset2,
                                                                begin_offset2 +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                begin_offset3,
                                                                begin_offset3 +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[0]],stride * dimension_offsets[dims[1]], stride * dimension_offsets[dims[2]],interp_func, pb,tuning);
            }*/




            return predict_error;
        }
        bool anchor=false;
        int interpolation_level = -1;
        uint blocksize;
        /*
        uint8_t interpolator_id;
        uint8_t interp_paradigm;
        uint8_t cubicSplineType=0;
        uint8_t direction_sequence_id;
        uint8_t adj_interp=0;
        */
        QoZ::Interp_Meta interp_meta;

        double eb_ratio = 0.5;
        double alpha;
        double beta;
        std::vector<std::string> interpolators = {"linear", "cubic"};
        std::vector<int> quant_inds;
        std::vector<bool> mark;
        size_t quant_index = 0; // for decompress
        size_t maxStep=0;
        //double max_error;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        size_t num_elements;

        std::array<size_t, N> global_dimensions;
        std::array<size_t, N> dimension_offsets;
        std::vector<std::array<int, N>> dimension_sequences;
        

        std::vector<float> prediction_errors;//for test, to delete in final version. The float time is to match the vector in config.
        //int peTracking=0;//for test, to delete in final version

        size_t cur_level; //temp for "adaptive anchor stride";
        //size_t min_anchor_level;//temp for "adaptive anchor stride";
       // double anchor_threshold=0.0;//temp for "adaptive anchor stride";



    };


};


#endif

