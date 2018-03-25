classdef TrackingErrors
    %TrackingErrors contains a set of enums to describe potential errors
    %in event-based feature tracking.
    %
    % Author: Alex Zihao Zhu, University of Pennsylvania
    % Email: alexzhu(at)seas.upenn.edu
    % Copyright 2018 University of Pennsylvania
    % Alex Zihao Zhu, Nikolay Atanasov, Kostas Daniilidis
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS, CONTRIBUTORS, AND THE
    % TRUSTEES OF THE UNIVERSITY OF PENNSYLVANIA "AS IS" AND ANY EXPRESS OR
    % IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
    % OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
    % IN NO EVENT SHALL THE COPYRIGHT OWNER, CONTRIBUTORS OR THE TRUSTEES OF
    % THE UNIVERSITY OF PENNSYLVANIA BE LIABLE FOR ANY DIRECT, INDIRECT,
    % INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
    % NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    % DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    % THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    % (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
    % THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    enumeration
        None, % No error.
        OutOfImg, % Feature left the image.
        Flow, % EM1_flow failed.
        ICPMaxIters, % EM2_affine took too many iterations to converge.
        ICPError, % The final cost in EM2_affine was too high.
        ICPFailed, % EM2_affine failed for another reason.
        NotValid % The feature was already not valid when asked to compute flow.
    end
    
end

